/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2022 Kévin Dietrich.
 * All rights reserved.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

#include "points_depuis_vdb.hh"

#include <openvdb/math/Math.h>
#include <openvdb/math/Stencils.h>
#include <openvdb/points/PointDelete.h>
#include <openvdb/points/PointGroup.h>
#include <openvdb/points/PointScatter.h>
#include <openvdb/tools/GridOperators.h>  // for tools::cpt()
#include <openvdb/tools/Interpolation.h>  // for tools::BoxSampler
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/Morphology.h>  // for tools::dilateActiveValues()
#include <openvdb/tools/PointScatter.h>
#include <openvdb/tree/LeafManager.h>

#include "biblinternes/moultfilage/boucle.hh"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"
#include "outils.hh"

namespace kvdb {

////////////////////////////////////////

// Simple wrapper class required by openvdb::tools::UniformPointScatter and
// NonUniformPointScatter
class PointAccessor {
  public:
    PointAccessor(AdaptriceMaillageVDB &sortie_points_) : sortie_points(sortie_points_)
    {
    }

    void add(const openvdb::Vec3R &pos)
    {
        sortie_points.ajouteUnPoint(pos.x(), pos.y(), pos.z());
    }

  protected:
    AdaptriceMaillageVDB &sortie_points;
};

////////////////////////////////////////

/// @brief Functor to translate points toward an isosurface
class SnapPointsOp {
  public:
    using Sampler = openvdb::tools::BoxSampler;

    enum class TypeDePoints { Invalide, Natif, VDB };

    // Constructor for Houdini points
    SnapPointsOp(AdaptriceMaillageVDB &detail,
                 const PlageDonnees &range,
                 float spread,
                 float isovalue,
                 bool rebuild,
                 bool dilate,
                 openvdb::BoolGrid::Ptr mask,
                 InterruptriceVDB *interrupter)
        : mTypeDePoints(est_valide(range) && !range.empty() ? TypeDePoints::Natif :
                                                              TypeDePoints::Invalide),
          mDetail(&detail), mRange(range), mSpread(spread), mIsovalue(isovalue), mRebuild(rebuild),
          mDilate(dilate), mMask(mask), mBoss(interrupter)
    {
    }

    // Constructor for VDB points
    SnapPointsOp(openvdb::points::PointDataGrid &vdbpts,
                 float spread,
                 float isovalue,
                 bool rebuild,
                 bool dilate,
                 openvdb::BoolGrid::Ptr mask,
                 InterruptriceVDB *interrupter)
        : mVdbPoints(&vdbpts), mSpread(spread), mIsovalue(isovalue), mRebuild(rebuild),
          mDilate(dilate), mMask(mask), mBoss(interrupter)
    {
        const auto leafIter = vdbpts.tree().cbeginLeaf();
        const auto descriptor = leafIter->attributeSet().descriptor();
        mAttrIdx = descriptor.find("P");
        mTypeDePoints = (mAttrIdx != openvdb::points::AttributeSet::INVALID_POS) ?
                            TypeDePoints::VDB :
                            TypeDePoints::Invalide;
    }

    template <typename GridT>
    void operator()(const GridT &aGrid)
    {
        if (mTypeDePoints == TypeDePoints::Invalide)
            return;

        const GridT *grid = &aGrid;

        // Replace the input grid with a rebuilt narrow-band level set, if requested
        // (typically because the isovalue is nonzero).
        typename GridT::Ptr sdf;
        if (mRebuild) {
            const float width = openvdb::LEVEL_SET_HALF_WIDTH;
            sdf = openvdb::tools::levelSetRebuild(*grid,
                                                  mIsovalue,
                                                  /*exterior=*/width,
                                                  /*interior=*/width,
                                                  /*xform=*/nullptr,
                                                  mBoss);
            if (sdf) {
                grid = sdf.get();
                mMask.reset();  // no need for a mask now that the input is a narrow-band level set
            }
        }

        // Compute the closest point transform of the SDF.
        const auto cpt = [&]() {
            if (!mMask) {
                return openvdb::tools::cpt(*grid, /*threaded=*/true, mBoss);
            }
            else {
                if (mDilate) {
                    // Dilate the isosurface mask to produce a suitably large CPT mask,
                    // to avoid unnecessary work in case the input is a dense SDF.
                    const int iterations = static_cast<int>(openvdb::LEVEL_SET_HALF_WIDTH);
                    openvdb::tools::dilateActiveValues(mMask->tree(),
                                                       iterations,
                                                       openvdb::tools::NN_FACE_EDGE,
                                                       openvdb::tools::IGNORE_TILES);
                }
                return openvdb::tools::cpt(*grid, *mMask, /*threaded=*/true, mBoss);
            }
        }();

        const auto &xform = aGrid.transform();
        if (mTypeDePoints == TypeDePoints::Natif) {
            // Translate Houdini points toward the isosurface.
            boucle_parallele_legere(mRange, [&](const PlageDonnees &r) {
                const auto cptAcc = cpt->getConstAccessor();
                auto start = r.begin(), end = r.end();

                if (mBoss && mBoss->wasInterrupted())
                    return;

                for (auto offset = start; offset < end; ++offset) {
                    openvdb::Vec3d p{mDetail->getPoint(offset)};
                    // Compute the closest surface point by linear interpolation.
                    const auto surfaceP = Sampler::sample(cptAcc, xform.worldToIndex(p));
                    // Translate the input point toward the surface.
                    p = surfaceP + mSpread * (p - surfaceP);  // (1-spread)*surfaceP + spread*p
                    mDetail->setPoint(offset, p);
                }
            });
        }
        else /*if (mTypeDePoints == TypeDePoints::VDB)*/ {
            // Translate VDB points toward the isosurface.
            using LeafMgr = openvdb::tree::LeafManager<openvdb::points::PointDataTree>;
            LeafMgr leafMgr(mVdbPoints->tree());
            boucle_parallele_legere(leafMgr.leafRange(), [&](const LeafMgr::LeafRange &range) {
                const auto cptAcc = cpt->getConstAccessor();
                for (auto leafIter = range.begin(); leafIter; ++leafIter) {
                    if (mBoss && mBoss->wasInterrupted())
                        break;
                    // Get a handle to this leaf node's point position array.
                    auto &posArray = leafIter->attributeArray(mAttrIdx);
                    openvdb::points::AttributeWriteHandle<openvdb::Vec3f> posHandle(posArray);
                    // For each point in this leaf node...
                    for (auto idxIter = leafIter->beginIndexOn(); idxIter; ++idxIter) {
                        // The point position is in index space and is relative to
                        // the center of the voxel.
                        const auto idxCenter = idxIter.getCoord().asVec3d();
                        const auto idxP = posHandle.get(*idxIter) + idxCenter;
                        // Compute the closest surface point by linear interpolation.
                        const openvdb::Vec3f surfaceP(Sampler::sample(cptAcc, idxP));
                        // Translate the input point toward the surface.
                        auto p = xform.indexToWorld(idxP);
                        p = surfaceP + mSpread * (p - surfaceP);  // (1-spread)*surfaceP + spread*p
                        // Transform back to index space relative to the voxel center.
                        posHandle.set(*idxIter, xform.worldToIndex(p) - idxCenter);
                    }
                }
            });
        }
    }

  private:
    TypeDePoints mTypeDePoints = TypeDePoints::Invalide;
    openvdb::points::PointDataGrid *mVdbPoints = nullptr;  // VDB points to be processed
    openvdb::Index64 mAttrIdx = openvdb::points::AttributeSet::INVALID_POS;
    AdaptriceMaillageVDB *mDetail =
        nullptr;                // the detail containing Houdini points to be processed
    PlageDonnees mRange{0, 0};  // the range of points to be processed
    float mSpread = 1;          // if 0, place points on the isosurface; if 1, don't move them
    float mIsovalue = 0;
    bool mRebuild = false;         // if true, generate a new SDF from the input grid
    bool mDilate = false;          // if true, dilate the isosurface mask
    openvdb::BoolGrid::Ptr mMask;  // an optional isosurface mask
    InterruptriceVDB *mBoss = nullptr;
};  // class SnapPointsOp

////////////////////////////////////////

struct BaseScatter {
    using NullCodec = openvdb::points::NullCodec;
    using FixedCodec16 = openvdb::points::FixedPointCodec<false>;
    using FixedCodec8 = openvdb::points::FixedPointCodec<true>;

    using PositionArray = openvdb::points::TypedAttributeArray<openvdb::Vec3f, NullCodec>;
    using PositionArray16 = openvdb::points::TypedAttributeArray<openvdb::Vec3f, FixedCodec16>;
    using PositionArray8 = openvdb::points::TypedAttributeArray<openvdb::Vec3f, FixedCodec8>;

    BaseScatter(const unsigned int seed, const float spread, InterruptriceVDB *interrupter)
        : mPoints(), mSeed(seed), mSpread(spread), mInterrupter(interrupter)
    {
    }
    virtual ~BaseScatter()
    {
    }

    /// @brief Print information about the scattered points
    /// @parm name  A name to insert into the printed info
    /// @parm os    The output stream
    virtual void print(const std::string &name, std::ostream &os = std::cout) const
    {
        if (!mPoints)
            return;
        const openvdb::Index64 points = openvdb::points::pointCount(mPoints->tree());
        const openvdb::Index64 voxels = mPoints->activeVoxelCount();
        os << points << " points into " << voxels << " active voxels in \"" << name
           << "\" corresponding to " << (double(points) / double(voxels)) << " points per voxel."
           << std::endl;
    }

    inline openvdb::points::PointDataGrid::Ptr points()
    {
        assert(mPoints);
        return mPoints;
    }

  protected:
    openvdb::points::PointDataGrid::Ptr mPoints;
    const unsigned int mSeed;
    const float mSpread;
    InterruptriceVDB *mInterrupter;
};  // BaseScatter

struct VDBUniformScatter : public BaseScatter {
    VDBUniformScatter(const openvdb::Index64 count,
                      const unsigned int seed,
                      const float spread,
                      const int compression,
                      InterruptriceVDB *interrupter)
        : BaseScatter(seed, spread, interrupter), mCount(count), mCompression(compression)
    {
    }

    template <typename PositionT, typename GridT>
    inline void resolveCompression(const GridT &grid)
    {
        using namespace openvdb::points;
        using PointDataGridT =
            openvdb::Grid<typename TreeConverter<typename GridT::TreeType>::Type>;
        mPoints =
            openvdb::points::uniformPointScatter<GridT, std::mt19937, PositionT, PointDataGridT>(
                grid, mCount, mSeed, mSpread, mInterrupter);
    }

    template <typename GridT>
    inline void operator()(const GridT &grid)
    {
        if (mCompression == 1) {
            this->resolveCompression<PositionArray16>(grid);
        }
        else if (mCompression == 2) {
            this->resolveCompression<PositionArray8>(grid);
        }
        else {
            this->resolveCompression<PositionArray>(grid);
        }
    }

    void print(const std::string &name, std::ostream &os = std::cout) const override
    {
        os << "Uniformly scattered ";
        BaseScatter::print(name, os);
    }

    const openvdb::Index64 mCount;
    const int mCompression;
};  // VDBUniformScatter

struct VDBDenseUniformScatter : public BaseScatter {
    VDBDenseUniformScatter(const float pointsPerVoxel,
                           const unsigned int seed,
                           const float spread,
                           const int compression,
                           InterruptriceVDB *interrupter)
        : BaseScatter(seed, spread, interrupter), mPointsPerVoxel(pointsPerVoxel),
          mCompression(compression)
    {
    }

    template <typename PositionT, typename GridT>
    inline void resolveCompression(const GridT &grid)
    {
        using namespace openvdb::points;
        using PointDataGridT =
            openvdb::Grid<typename TreeConverter<typename GridT::TreeType>::Type>;
        mPoints = denseUniformPointScatter<GridT, std::mt19937, PositionT, PointDataGridT>(
            grid, mPointsPerVoxel, mSeed, mSpread, mInterrupter);
    }

    template <typename GridT>
    inline void operator()(const GridT &grid)
    {
        if (mCompression == 1) {
            this->resolveCompression<PositionArray16>(grid);
        }
        else if (mCompression == 2) {
            this->resolveCompression<PositionArray8>(grid);
        }
        else {
            this->resolveCompression<PositionArray>(grid);
        }
    }

    void print(const std::string &name, std::ostream &os = std::cout) const override
    {
        os << "Dense uniformly scattered ";
        BaseScatter::print(name, os);
    }

    const float mPointsPerVoxel;
    const int mCompression;
};  // VDBDenseUniformScatter

struct VDBNonUniformScatter : public BaseScatter {
    VDBNonUniformScatter(const float pointsPerVoxel,
                         const unsigned int seed,
                         const float spread,
                         const int compression,
                         InterruptriceVDB *interrupter)
        : BaseScatter(seed, spread, interrupter), mPointsPerVoxel(pointsPerVoxel),
          mCompression(compression)
    {
    }

    template <typename PositionT, typename GridT>
    inline void resolveCompression(const GridT &grid)
    {
        using namespace openvdb::points;
        using PointDataGridT =
            openvdb::Grid<typename TreeConverter<typename GridT::TreeType>::Type>;
        mPoints = nonUniformPointScatter<GridT, std::mt19937, PositionT, PointDataGridT>(
            grid, mPointsPerVoxel, mSeed, mSpread, mInterrupter);
    }

    template <typename GridT>
    inline void operator()(const GridT &grid)
    {
        if (mCompression == 1) {
            this->resolveCompression<PositionArray16>(grid);
        }
        else if (mCompression == 2) {
            this->resolveCompression<PositionArray8>(grid);
        }
        else {
            this->resolveCompression<PositionArray>(grid);
        }
    }

    void print(const std::string &name, std::ostream &os = std::cout) const override
    {
        os << "Non-uniformly scattered ";
        BaseScatter::print(name, os);
    }

    const float mPointsPerVoxel;
    const int mCompression;
};  // VDBNonUniformScatter

template <typename SurfaceGridT>
struct MarkPointsOutsideIso {
    using GroupIndex = openvdb::points::AttributeSet::Descriptor::GroupIndex;
    using LeafManagerT = openvdb::tree::LeafManager<openvdb::points::PointDataTree>;
    using PositionHandleT =
        openvdb::points::AttributeHandle<openvdb::Vec3f, openvdb::points::NullCodec>;
    using SurfaceValueT = typename SurfaceGridT::ValueType;

    MarkPointsOutsideIso(const SurfaceGridT &grid, const GroupIndex &deadIndex)
        : mGrid(grid), mDeadIndex(deadIndex)
    {
    }

    void operator()(const LeafManagerT::LeafRange &range) const
    {
        openvdb::math::BoxStencil<const SurfaceGridT> stencil(mGrid);
        for (auto leaf = range.begin(); leaf; ++leaf) {

            PositionHandleT::Ptr positionHandle = PositionHandleT::create(
                leaf->constAttributeArray(0));
            openvdb::points::GroupWriteHandle deadHandle = leaf->groupWriteHandle(mDeadIndex);

            for (auto voxel = leaf->cbeginValueOn(); voxel; ++voxel) {

                const openvdb::Coord &ijk = voxel.getCoord();
                const openvdb::Vec3d vec = ijk.asVec3d();

                for (auto iter = leaf->beginIndexVoxel(ijk); iter; ++iter) {
                    const openvdb::Index index = *iter;
                    const openvdb::Vec3d pos = openvdb::Vec3d(positionHandle->get(index)) + vec;

                    stencil.moveTo(pos);
                    if (stencil.interpolation(pos) > openvdb::zeroVal<SurfaceValueT>()) {
                        deadHandle.set(index, true);
                    }
                }
            }
        }
    }

  private:
    const SurfaceGridT &mGrid;
    const GroupIndex &mDeadIndex;
};  // MarkPointsOutsideIso

template <typename OpType>
inline bool process(const openvdb::GridBase &grid, OpType &op, const std::string *name)
{
    bool success = grid.apply<AllGridTypes>(op);
    if (name)
        op.print(*name);
    return success;
}

// Extract an SDF interior mask in which to scatter points.
inline openvdb::BoolGrid::Ptr extractInteriorMask(const openvdb::GridBase::ConstPtr grid,
                                                  const float isovalue)
{
    if (grid->isType<openvdb::FloatGrid>()) {
        return openvdb::tools::sdfInteriorMask(static_cast<const openvdb::FloatGrid &>(*grid),
                                               isovalue);
    }
    else if (grid->isType<openvdb::DoubleGrid>()) {
        return openvdb::tools::sdfInteriorMask(static_cast<const openvdb::DoubleGrid &>(*grid),
                                               isovalue);
    }
    return nullptr;
}

// Extract an SDF isosurface mask in which to scatter points.
inline openvdb::BoolGrid::Ptr extractIsosurfaceMask(const openvdb::GridBase::ConstPtr grid,
                                                    const float isovalue)
{
    if (grid->isType<openvdb::FloatGrid>()) {
        return openvdb::tools::extractIsosurfaceMask(
            static_cast<const openvdb::FloatGrid &>(*grid), isovalue);
    }
    else if (grid->isType<openvdb::DoubleGrid>()) {
        return openvdb::tools::extractIsosurfaceMask(
            static_cast<const openvdb::DoubleGrid &>(*grid), double(isovalue));
    }
    return nullptr;
}

// Remove VDB Points scattered outside of a level set
inline void cullVDBPoints(openvdb::points::PointDataTree &tree,
                          const openvdb::GridBase::ConstPtr grid)
{
    const auto leaf = tree.cbeginLeaf();
    if (leaf) {
        using GroupIndex = openvdb::points::AttributeSet::Descriptor::GroupIndex;
        openvdb::points::appendGroup(tree, "dead");
        const GroupIndex idx = leaf->attributeSet().groupIndex("dead");

        openvdb::tree::LeafManager<openvdb::points::PointDataTree> leafManager(tree);

        if (grid->isType<openvdb::FloatGrid>()) {
            const openvdb::FloatGrid &typedGrid = static_cast<const openvdb::FloatGrid &>(*grid);
            MarkPointsOutsideIso<openvdb::FloatGrid> mark(typedGrid, idx);
            tbb::parallel_for(leafManager.leafRange(), mark);
        }
        else if (grid->isType<openvdb::DoubleGrid>()) {
            const openvdb::DoubleGrid &typedGrid = static_cast<const openvdb::DoubleGrid &>(*grid);
            MarkPointsOutsideIso<openvdb::DoubleGrid> mark(typedGrid, idx);
            tbb::parallel_for(leafManager.leafRange(), mark);
        }
        openvdb::points::deleteFromGroup(tree, "dead");
    }
}

void points_depuis_vdb(ContexteKuri &ctx_kuri,
                       EnveloppeContexteEvaluation &ctx_eval,
                       ParametresPointsDepuisVDB &params,
                       InterruptriceVDB &boss)
{
    boss.start("Création de points dans des grilles OpenVDB");

    if (!params.groupe_grilles) {
        ctx_eval.rapporteErreur("Aucun groupe de grilles en entrée.");
        return;
    }

    const int seed = params.graine;
    const auto theSpread = params.diffusion;
    const bool verbose = params.verbeux;
    const openvdb::Index64 pointCount = params.nombre_de_points;
    const float ptsPerVox = params.points_par_voxel;
    const auto sdfdomain = params.domaine_champs_de_distance;
    const float density = params.densite;
    const bool multiplyDensity = params.multiplie_par_densite;
    const auto theIsovalue = params.valeur_iso;
    const PolitiqueNommageGrille outputName = params.nommage_grilles_sortie;
    const std::string customName = outils::chaine_depuis_accesseuse(params.nom_grille_sortie);

    // Get the group of grids to process.
    const auto group = outils::grilles_depuis_iteratrice(*params.groupe_grilles);
    if (group.empty()) {
        ctx_eval.rapporteErreur("Aucune grille trouvée dans le groupe de grilles.");
        return;
    }

    // Choose a fast random generator with a long period. Drawback here for
    // mt11213b is that it requires 352*sizeof(uint32) bytes.
    using RandGen = std::mersenne_twister_engine<uint32_t,
                                                 32,
                                                 351,
                                                 175,
                                                 19,
                                                 0xccab8ee7,
                                                 11,
                                                 0xffffffff,
                                                 7,
                                                 0x31b6ab00,
                                                 15,
                                                 0xffe50000,
                                                 17,
                                                 1812433253>;  // mt11213b
    RandGen mtRand(seed);

    const auto pmode = params.mode_generation_points;
    const bool vdbPoints = params.genere_grille_points;
    const bool clipPoints = vdbPoints && params.rogne_selon_surface;
    const int posCompression = vdbPoints ? params.compression_desiree : 0;
    const bool snapPointsToSurface = ((sdfdomain == SURFACE) &&
                                      !openvdb::math::isApproxEqual(theSpread, 1.0f));

    // If the domain is the isosurface, set the spread to 1 while generating points
    // so that each point ends up snapping to a unique point on the surface.
    const float spread = (snapPointsToSurface ? 1.f : theSpread);

    std::vector<std::string> emptyGrids;
    std::vector<openvdb::points::PointDataGrid::Ptr> pointGrids;
    AdaptriceMaillageVDB sortie_points = AdaptriceMaillageVDB::enveloppe(params.sortie_points);
    PointAccessor pointAccessor(sortie_points);

    // Process each VDB primitive (with a non-null grid pointer)
    // that belongs to the selected group.
    for (auto primIter : group) {
        // Retrieve a read-only grid pointer.
        openvdb::GridBase::ConstPtr grid = primIter->grid;
        TypeVolume gridType = VDB_type_volume_pour_grille(primIter);
        const std::string gridName = grid->getName();

        if (grid->empty()) {
            emptyGrids.push_back(gridName);
            continue;
        }

        const std::string *const name = verbose ? &gridName : nullptr;
        const openvdb::GridClass gridClass = grid->getGridClass();
        const bool isSignedDistance = (gridClass == openvdb::GRID_LEVEL_SET);
        bool performCull = false;

        const auto isovalue = (gridClass != openvdb::GRID_FOG_VOLUME) ?
                                  theIsovalue :
                                  openvdb::math::Clamp(
                                      theIsovalue, openvdb::math::Tolerance<float>::value(), 1.f);

        openvdb::BoolGrid::Ptr mask;
        if (sdfdomain != BANDE_FINE) {
            auto iso = isovalue;
            if (clipPoints) {
                const openvdb::Vec3d voxelSize = grid->voxelSize();
                const double maxVoxelSize = openvdb::math::Max(
                    voxelSize.x(), voxelSize.y(), voxelSize.z());
                iso += static_cast<float>(maxVoxelSize / 2.0);
                performCull = true;
            }

            if (sdfdomain == INTERIEUR) {
                if (isSignedDistance) {
                    // If the input is an SDF, compute a mask of its interior.
                    // (Fog volumes are their own interior masks.)
                    mask = extractInteriorMask(grid, iso);
                }
            }
            else if (sdfdomain == SURFACE) {
                mask = extractIsosurfaceMask(grid, iso);
            }
            if (mask) {
                grid = mask;
                gridType = VOLUME_BOOL;
            }
        }

        std::string vdbName;
        if (vdbPoints) {
            if (outputName == GARDE_NOM_ORIGINAL)
                vdbName = gridName;
            else if (outputName == AJOUTE_SUFFIXE)
                vdbName = gridName + customName;
            else
                vdbName = customName;
        }

        openvdb::points::PointDataGrid::Ptr pointGrid;

        const auto postprocessVDBPoints = [&](BaseScatter &scatter, bool cull) {
            pointGrid = scatter.points();
            if (cull) {
                cullVDBPoints(pointGrid->tree(), grid);
            }
            pointGrid->setName(vdbName);
            pointGrids.push_back(pointGrid);
            if (verbose)
                scatter.print(gridName);
        };

        using DenseScatterer =
            openvdb::tools::DenseUniformPointScatter<PointAccessor, RandGen, InterruptriceVDB>;
        using NonuniformScatterer =
            openvdb::tools::NonUniformPointScatter<PointAccessor, RandGen, InterruptriceVDB>;
        using UniformScatterer =
            openvdb::tools::UniformPointScatter<PointAccessor, RandGen, InterruptriceVDB>;

        const auto startOffset = sortie_points.pointCount();

        switch (pmode) {

            case PAR_COMPTE_TOTAL:
                if (vdbPoints) {  // VDB points
                    VDBUniformScatter scatter(pointCount, seed, spread, posCompression, &boss);
                    if (process(*grid, scatter, name)) {
                        postprocessVDBPoints(scatter, performCull);
                    }
                }
                else {  // Houdini points
                    UniformScatterer scatter(pointAccessor, pointCount, mtRand, spread, &boss);
                    process(*grid, scatter, name);
                }
                break;

            case PAR_DENSITE_LOCALE:
                if (multiplyDensity && !isSignedDistance) {  // local density
                    if (vdbPoints) {                         // VDB points
                        const auto dim = grid->transform().voxelSize();
                        VDBNonUniformScatter scatter(static_cast<float>(density * dim.product()),
                                                     seed,
                                                     spread,
                                                     posCompression,
                                                     &boss);
                        if (!grid->apply<openvdb::NumericGridTypes>(scatter)) {
                            throw std::runtime_error(
                                "Only scalar grids support voxel scaling of density");
                        }
                        postprocessVDBPoints(scatter, /*cull=*/false);
                    }
                    else {  // Houdini points
                        NonuniformScatterer scatter(pointAccessor, density, mtRand, spread, &boss);
                        if (!grid->apply<openvdb::NumericGridTypes>(scatter)) {
                            throw std::runtime_error(
                                "Only scalar grids support voxel scaling of density");
                        }
                        if (verbose)
                            scatter.print(gridName);
                    }
                }
                else {                // global density
                    if (vdbPoints) {  // VDB points
                        const auto dim = grid->transform().voxelSize();
                        const auto totalPointCount = openvdb::Index64(
                            density * dim.product() * double(grid->activeVoxelCount()));
                        VDBUniformScatter scatter(
                            totalPointCount, seed, spread, posCompression, &boss);
                        if (process(*grid, scatter, name)) {
                            postprocessVDBPoints(scatter, performCull);
                        }
                    }
                    else {  // Houdini points
                        UniformScatterer scatter(pointAccessor, density, mtRand, spread, &boss);
                        process(*grid, scatter, name);
                    }
                }
                break;

            case PAR_VOXEL:
                if (vdbPoints) {  // VDB points
                    VDBDenseUniformScatter scatter(ptsPerVox, seed, spread, posCompression, &boss);
                    if (process(*grid, scatter, name)) {
                        postprocessVDBPoints(scatter, performCull);
                    }
                }
                else {  // Houdini points
                    DenseScatterer scatter(pointAccessor, ptsPerVox, mtRand, spread, &boss);
                    process(*grid, scatter, name);
                }
                break;

            default:
                throw std::runtime_error("Expected 0, 1 or 2 for \"pointmode\", got " +
                                         std::to_string(pmode));
        }  // switch pmode

        if (snapPointsToSurface) {
            // Dilate the mask if it is a single-voxel-wide isosurface mask.
            const bool dilate = (mask && (sdfdomain == SURFACE));
            // Generate a new SDF if the input is a fog volume or if the isovalue is nonzero.
            const bool rebuild = (!isSignedDistance || !openvdb::math::isApproxZero(isovalue));
            if (!vdbPoints) {
                const auto range = sortie_points.plagePoint(startOffset);
                // Use the original spread value to control how close to the surface points lie.
                SnapPointsOp op{
                    sortie_points, range, theSpread, isovalue, rebuild, dilate, mask, &boss};
                grid->apply<openvdb::RealGridTypes>(op);
            }
            else if (vdbPoints && pointGrid) {
                SnapPointsOp op{*pointGrid, theSpread, isovalue, rebuild, dilate, mask, &boss};
                grid->apply<openvdb::RealGridTypes>(op);
            }
        }
    }  // for each grid

    if (!emptyGrids.empty()) {
        std::string s = "Les grilles suivantes étaient vides : " +
                        outils::enchaine(emptyGrids, ", ");
        ctx_eval.rapporteAvertissement(s);
    }

    /* Crée un groupe pour les points au besoin. */
    if (params.cree_groupe) {
        const std::string groupName = outils::chaine_depuis_accesseuse(params.nom_groupe_sortie);
        auto groupe_de_points = sortie_points.creeUnGroupeDePoints(groupName);

        if (groupe_de_points) {
            sortie_points.ajoutePlageAuGroupe(groupe_de_points, 0, sortie_points.pointCount());
        }

        for (auto &pointGrid : pointGrids) {
            openvdb::points::appendGroup(pointGrid->tree(), groupName);
            openvdb::points::setGroup(pointGrid->tree(), groupName);
        }
    }

    for (auto &pointGrid : pointGrids) {
        outils::exporte_grille_vdb(
            &ctx_kuri, params.sortie_grilles, pointGrid, pointGrid->getName().c_str());
    }

    boss.end();
}

}  // namespace kvdb
