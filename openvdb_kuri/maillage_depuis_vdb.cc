/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "maillage_depuis_vdb.hh"

#include <list>

#include "openvdb/math/Ray.h"
#include "openvdb/tools/Mask.h"  // pour tools::interiorMask()
#include "openvdb/tools/Morphology.h"
#include "openvdb/tools/VolumeToMesh.h"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"
#include "outils.hh"
#include "vdb_depuis_maillage.hh"

using namespace openvdb;

namespace kvdb {

struct InteriorMaskOp {
    InteriorMaskOp(double iso = 0.0) : inIsovalue(iso)
    {
    }

    template <typename GridType>
    void operator()(const GridType &grid)
    {
        outGridPtr = tools::interiorMask(grid, inIsovalue);
    }

    const double inIsovalue;
    BoolGrid::Ptr outGridPtr;
};

// Extract a boolean mask from a grid of any type.
inline GridBase::ConstPtr getMaskFromGrid(const GridBase::ConstPtr &gridPtr, double isovalue = 0.0)
{
    if (!gridPtr) {
        return nullptr;
    }

    if (gridPtr->isType<BoolGrid>()) {
        // If the input grid is already boolean, return it.
        return gridPtr;
    }

    InteriorMaskOp op{isovalue};
    gridPtr->apply<AllGridTypes>(op);
    return op.outGridPtr;
}

static void ajoute_masque_surface(EnveloppeContexteEvaluation &ctx_eval,
                                  ParametresVDBVersMaillage *params,
                                  tools::VolumeToMesh &mesher)
{
    if (!params->grille_masque_surface) {
        return;
    }

    if (!params->grille_masque_surface->grid) {
        ctx_eval.rapporteAvertissement("Aucune grille pour le masque de surface");
        return;
    }

    auto grille_masque = getMaskFromGrid(params->grille_masque_surface->grid,
                                         params->decalage_masque);

    mesher.setSurfaceMask(grille_masque, params->inverse_masque);
}

static void ajoute_champs_adaptivite(EnveloppeContexteEvaluation &ctx_eval,
                                     ParametresVDBVersMaillage *params,
                                     tools::VolumeToMesh &mesher)
{
    if (!params->grille_champs_adaptivite) {
        return;
    }

    if (!params->grille_champs_adaptivite->grid) {
        ctx_eval.rapporteAvertissement("Aucune grille pour le champs d'adaptivité");
        return;
    }

    if (VDB_type_volume_pour_grille(params->grille_champs_adaptivite) != VOLUME_R32) {
        ctx_eval.rapporteAvertissement("La grille du champs d'adaptivité n'est pas de type réel");
        return;
    }

    auto grille = gridConstPtrCast<FloatGrid>(params->grille_champs_adaptivite->grid);
    mesher.setSpatialAdaptivity(grille);
}

enum {
    POLYGONE_INTERNE = 0,
    POLYGONE_INTERNE_SUR_COUTURE = 1,
    POLYGONE_SUPERFICIEL = 2,
    POLYGONE_SUPERFICIEL_SUR_COUTURE = 3,

    NOMBRE_TYPE_POLYGONE = 4,
};

/* Convertis les drapeaux du polygone en une valeur de l'énumération ci-dessus. */
static inline int convertis_drapeaux_polygone(char drapeau_vdb)
{
    constexpr char est_externe = char(tools::POLYFLAG_EXTERIOR);
    constexpr char est_sur_couture = char(tools::POLYFLAG_FRACTURE_SEAM);
    return (((drapeau_vdb & est_externe) != 0) << 1) | ((drapeau_vdb & est_sur_couture) != 0);
}

void copyMesh(AdaptriceMaillageVDB &maillage, tools::VolumeToMesh &mesher, const char *gridName)
{
    /* Exporte les points. */
    const tools::PointList &points = mesher.pointList();

    const int decalage_point = maillage.pointCount();
    maillage.ajoutePoints(reinterpret_cast<float *>(points.get()),
                          static_cast<long>(mesher.pointListSize()));

    if (mesher.pointFlags().size() == mesher.pointListSize()) {
        void *groupe_points_sur_couture = maillage.creeUnGroupeDePoints("points_couture");
        if (groupe_points_sur_couture) {
            for (int i = 0; i < mesher.pointListSize(); i++) {
                if (mesher.pointFlags()[i]) {
                    maillage.ajouteAuGroupe(groupe_points_sur_couture, i);
                }
            }
        }
    }

    /* Exporte les polygones. */
    tools::PolygonPoolList &polygonPoolList = mesher.polygonPoolList();

    long nquads[NOMBRE_TYPE_POLYGONE] = {0, 0, 0, 0};
    long ntris[NOMBRE_TYPE_POLYGONE] = {0, 0, 0, 0};
    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
        const tools::PolygonPool &polygons = polygonPoolList[n];
        for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
            int flags = convertis_drapeaux_polygone(polygons.quadFlags(i));
            ++nquads[flags];
        }
        for (size_t i = 0, I = polygons.numTriangles(); i < I; ++i) {
            int flags = convertis_drapeaux_polygone(polygons.triangleFlags(i));
            ++ntris[flags];
        }
    }

    long nverts[NOMBRE_TYPE_POLYGONE] = {nquads[0] * 4 + ntris[0] * 3,
                                         nquads[1] * 4 + ntris[1] * 3,
                                         nquads[2] * 4 + ntris[2] * 3,
                                         nquads[3] * 4 + ntris[3] * 3};
    std::vector<int> verts[NOMBRE_TYPE_POLYGONE];
    for (int flags = 0; flags < 4; ++flags) {
        verts[flags].resize(nverts[flags]);
    }

    long iquad[NOMBRE_TYPE_POLYGONE] = {0, 0, 0, 0};
    long itri[NOMBRE_TYPE_POLYGONE] = {nquads[0] * 4, nquads[1] * 4, nquads[2] * 4, nquads[3] * 4};

    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
        const tools::PolygonPool &polygons = polygonPoolList[n];

        // Copy quads
        for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
            const Vec4I &quad = polygons.quad(i);
            int flags = convertis_drapeaux_polygone(polygons.quadFlags(i));
            verts[flags][iquad[flags]++] = quad[0] + decalage_point;
            verts[flags][iquad[flags]++] = quad[1] + decalage_point;
            verts[flags][iquad[flags]++] = quad[2] + decalage_point;
            verts[flags][iquad[flags]++] = quad[3] + decalage_point;
        }

        // Copy triangles (adaptive mesh)
        for (size_t i = 0, I = polygons.numTriangles(); i < I; ++i) {
            const Vec3I &triangle = polygons.triangle(i);
            int flags = convertis_drapeaux_polygone(polygons.triangleFlags(i));
            verts[flags][itri[flags]++] = triangle[0] + decalage_point;
            verts[flags][itri[flags]++] = triangle[1] + decalage_point;
            verts[flags][itri[flags]++] = triangle[2] + decalage_point;
        }
    }

    void *groupe_polygones_sur_couture = maillage.creeUnGroupeDePolygones("polygones_sur_couture");
    void *groupe_polygones_sur_surface = maillage.creeUnGroupeDePolygones("polygones_sur_surface");
    void *groupe_polygones_internes = maillage.creeUnGroupeDePolygones("polygones_internes");

    std::vector<int> sommets_par_polygone;
    long decalage_groupe = maillage.polygonCount();
    for (int flags = 0; flags < 4; ++flags) {

        if (!nquads[flags] && !ntris[flags])
            continue;

        sommets_par_polygone.resize(nquads[flags] + ntris[flags]);
        std::fill(sommets_par_polygone.begin(), sommets_par_polygone.begin() + nquads[flags], 4);
        std::fill(sommets_par_polygone.begin() + nquads[flags],
                  sommets_par_polygone.begin() + nquads[flags] + ntris[flags],
                  3);

        maillage.ajouteListePolygones(verts[flags].data(),
                                      sommets_par_polygone.data(),
                                      static_cast<long>(sommets_par_polygone.size()));
        long fin_groupe = decalage_groupe + nquads[flags] + ntris[flags];

        if (groupe_polygones_sur_couture && (flags & 1)) {
            maillage.ajoutePlageAuGroupe(
                groupe_polygones_sur_couture, decalage_groupe, fin_groupe);
        }
        if (groupe_polygones_sur_surface && (flags & 2)) {
            maillage.ajoutePlageAuGroupe(
                groupe_polygones_sur_surface, decalage_groupe, fin_groupe);
        }
        if (groupe_polygones_internes && !(flags & 2)) {
            maillage.ajoutePlageAuGroupe(groupe_polygones_internes, decalage_groupe, fin_groupe);
        }

        decalage_groupe += nquads[flags] + ntris[flags];
    }

    // À FAIRE : attribut chaine Keep VDB grid name
    //    const GA_Index lastPrim = detail.getNumPrimitives();
    //    if (gridName != nullptr && firstPrim != lastPrim) {

    //        GA_RWAttributeRef aRef = detail.findPrimitiveAttribute("name");
    //        if (!aRef.isValid())
    //            aRef = detail.addStringTuple(GA_ATTRIB_PRIMITIVE, "name", 1);

    //        GA_Attribute *nameAttr = aRef.getAttribute();
    //        if (nameAttr) {
    //            const GA_AIFSharedStringTuple *stringAIF = nameAttr->getAIFSharedStringTuple();
    //            if (stringAIF) {
    //                GA_Range range(detail.getPrimitiveMap(),
    //                               detail.primitiveOffset(firstPrim),
    //                               detail.primitiveOffset(lastPrim));
    //                stringAIF->setString(nameAttr, range, gridName, 0);
    //            }
    //        }
    //    }
}

////////////////////////////////////////

/// TBB body object for threaded sharp feature construction
template <typename IndexTreeType, typename BoolTreeType>
class GenAdaptivityMaskOp {
  public:
    using BoolLeafManager = tree::LeafManager<BoolTreeType>;

    GenAdaptivityMaskOp(const AdaptriceMaillageVDB &refGeo,
                        const IndexTreeType &indexTree,
                        BoolLeafManager &,
                        float edgetolerance = 0.0);

    void run(bool threaded = true);

    void operator()(const tbb::blocked_range<size_t> &) const;

  private:
    const AdaptriceMaillageVDB &mRefGeo;
    const IndexTreeType &mIndexTree;
    BoolLeafManager &mLeafs;
    float mEdgeTolerance;
    std::vector<Vec3s> normaux;
};

template <typename IndexTreeType, typename BoolTreeType>
GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::GenAdaptivityMaskOp(
    const AdaptriceMaillageVDB &refGeo,
    const IndexTreeType &indexTree,
    BoolLeafManager &leafMgr,
    float edgetolerance)
    : mRefGeo(refGeo), mIndexTree(indexTree), mLeafs(leafMgr), mEdgeTolerance(edgetolerance)
{
    mEdgeTolerance = std::max(0.0f, mEdgeTolerance);
    mEdgeTolerance = std::min(1.0f, mEdgeTolerance);

    normaux.resize(refGeo.polygonCount());

    for (size_t i = 0; i < refGeo.polygonCount(); i++) {
        normaux[i] = refGeo.normalPolygone(i);
    }
}

template <typename IndexTreeType, typename BoolTreeType>
void GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::run(bool threaded)
{
    if (threaded) {
        tbb::parallel_for(mLeafs.getRange(), *this);
    }
    else {
        (*this)(mLeafs.getRange());
    }
}

template <typename IndexTreeType, typename BoolTreeType>
void GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::operator()(
    const tbb::blocked_range<size_t> &range) const
{
    using IndexAccessorType = typename tree::ValueAccessor<const IndexTreeType>;
    IndexAccessorType idxAcc(mIndexTree);

    Vec3s tmpN, normal;
    int tmpIdx;

    Coord ijk, nijk;
    typename BoolTreeType::LeafNodeType::ValueOnIter iter;

    for (size_t n = range.begin(); n < range.end(); ++n) {
        iter = mLeafs.leaf(n).beginValueOn();
        for (; iter; ++iter) {
            ijk = iter.getCoord();

            bool edgeVoxel = false;

            int idx = idxAcc.getValue(ijk);

            normal = normaux[idx];

            for (size_t i = 0; i < 18; ++i) {
                nijk = ijk + util::COORD_OFFSETS[i];
                if (idxAcc.probeValue(nijk, tmpIdx) && tmpIdx != idx) {
                    tmpN = normaux[tmpIdx];

                    if (normal.dot(tmpN) < mEdgeTolerance) {
                        edgeVoxel = true;
                        break;
                    }
                }
            }

            if (!edgeVoxel)
                iter.setValueOff();
        }
    }
}

////////////////////////////////////////

/// TBB body object for threaded sharp feature construction
class SharpenFeaturesOp {
  public:
    using EdgeData = tools::MeshToVoxelEdgeData;

    SharpenFeaturesOp(AdaptriceMaillageVDB &meshGeo,
                      const AdaptriceMaillageVDB &refGeo,
                      EdgeData &edgeData,
                      const math::Transform &xform,
                      const void *surfacePrims = nullptr,
                      const BoolTree *mask = nullptr);

    void operator()(const tbb::blocked_range<long> &) const;

  private:
    AdaptriceMaillageVDB &mMeshGeo;
    const AdaptriceMaillageVDB &mRefGeo;
    EdgeData &mEdgeData;
    const math::Transform &mXForm;
    const void *mSurfacePrims;
    const BoolTree *mMaskTree;
};

SharpenFeaturesOp::SharpenFeaturesOp(AdaptriceMaillageVDB &meshGeo,
                                     const AdaptriceMaillageVDB &refGeo,
                                     EdgeData &edgeData,
                                     const math::Transform &xform,
                                     const void *surfacePrims,
                                     const BoolTree *mask)
    : mMeshGeo(meshGeo), mRefGeo(refGeo), mEdgeData(edgeData), mXForm(xform),
      mSurfacePrims(surfacePrims), mMaskTree(mask)
{
}

void SharpenFeaturesOp::operator()(const tbb::blocked_range<long> &range) const
{
    tools::MeshToVoxelEdgeData::Accessor acc = mEdgeData.getAccessor();

    using BoolAccessor = tree::ValueAccessor<const BoolTree>;
    std::unique_ptr<BoolAccessor> maskAcc;

    if (mMaskTree) {
        maskAcc.reset(new BoolAccessor(*mMaskTree));
    }

    Vec3s tmpN, tmpP, avgP;
    math::BBox<Vec3d> cell;

    Vec3d pos, normal;
    Coord ijk;

    std::vector<Vec3d> points, normals;
    std::vector<Index32> primitives;

    points.reserve(12);
    normals.reserve(12);
    primitives.reserve(12);
    for (long ptnOffset = range.begin(); ptnOffset < range.end(); ++ptnOffset) {
        // Check if this point is referenced by a surface primitive.
        if (mSurfacePrims && !mMeshGeo.groupePolygonePossedePoint(mSurfacePrims, ptnOffset))
            continue;

        tmpP = mMeshGeo.getPoint(ptnOffset);
        pos[0] = tmpP.x();
        pos[1] = tmpP.y();
        pos[2] = tmpP.z();

        pos = mXForm.worldToIndex(pos);

        ijk[0] = int(std::floor(pos[0]));
        ijk[1] = int(std::floor(pos[1]));
        ijk[2] = int(std::floor(pos[2]));

        if (maskAcc && !maskAcc->isValueOn(ijk))
            continue;

        points.clear();
        normals.clear();
        primitives.clear();

        // get voxel-edge intersections
        mEdgeData.getEdgeData(acc, ijk, points, primitives);

        avgP = Vec3s(0.0f, 0.0f, 0.0f);

        // get normal list
        for (size_t n = 0, N = points.size(); n < N; ++n) {
            avgP.x() = static_cast<float>(avgP.x() + points[n].x());
            avgP.y() = static_cast<float>(avgP.y() + points[n].y());
            avgP.z() = static_cast<float>(avgP.z() + points[n].z());

            tmpN = mRefGeo.normalPolygone(primitives[n]);

            normal[0] = tmpN.x();
            normal[1] = tmpN.y();
            normal[2] = tmpN.z();

            normals.push_back(normal);
        }

        // Calculate feature point position
        if (points.size() <= 1) {
            continue;
        }

        pos = tools::findFeaturePoint(points, normals);

        // Constrain points to stay inside their initial
        // coordinate cell.
        auto min_bound = Vec3d(double(ijk[0]), double(ijk[1]), double(ijk[2]));
        auto max_bound = Vec3d(double(ijk[0] + 1), double(ijk[1] + 1), double(ijk[2] + 1));
        cell = math::BBox<Vec3d>(min_bound, max_bound);
        cell.expand(0.3);

        if (!cell.isInside(pos)) {
            Vec3s org(static_cast<float>(pos[0]),
                      static_cast<float>(pos[1]),
                      static_cast<float>(pos[2]));

            avgP *= 1.f / float(points.size());
            Vec3s dir = avgP - org;
            dir.normalize();

            double distance;

            math::Ray ray;
            ray.reset(org, dir);

            double t1;
            if (ray.intersects(cell, distance, t1)) {
                tmpP = org + dir * distance;

                pos[0] = tmpP.x();
                pos[1] = tmpP.y();
                pos[2] = tmpP.z();
            }
        }

        pos = mXForm.indexToWorld(pos);

        tmpP.x() = static_cast<float>(pos[0]);
        tmpP.y() = static_cast<float>(pos[1]);
        tmpP.z() = static_cast<float>(pos[2]);

        mMeshGeo.setPoint(ptnOffset, tmpP);
    }
}

template <typename GridType>
static void vdb_vers_polygones_reference(ContexteKuri *ctx,
                                         EnveloppeContexteEvaluation &ctx_eval,
                                         ParametresVDBVersMaillage *params,
                                         AdaptriceMaillageVDB &maillage,
                                         tools::VolumeToMesh &mesher,
                                         std::list<GridBase::ConstPtr> const &grids,
                                         InterruptriceVDB &boss)
{
    typename GridType::ConstPtr firstGrid = gridConstPtrCast<GridType>(grids.front());

    if (!firstGrid) {
        ctx_eval.rapporteErreur("Type de grille non supporté");
        return;
    }

    using TreeType = typename GridType::TreeType;
    using ValueType = typename GridType::ValueType;
    using IntGridT = typename GridType::template ValueConverter<Int32>::Type;
    math::Transform::Ptr transform = firstGrid->transform().copy();
    const ValueType backgroundValue = firstGrid->background();
    const GridClass gridClass = firstGrid->getGridClass();

    /* Crée la grille de référence. */
    tools::MeshToVoxelEdgeData edgeData;

    ExtractionDonneesVDBDepuisPolygones extraction;
    if (params->affiner_les_traits) {
        extraction.donnees_aretes_voxel = &edgeData;
    }

    float largeur_de_bande = 3.0f;
    if (gridClass != GRID_LEVEL_SET) {
        largeur_de_bande = static_cast<float>(backgroundValue) /
                           static_cast<float>(transform->voxelSize()[0]);
    }

    GrilleVDB grille_reference_maillage;
    grille_reference_maillage.grid = ConstPtrCast<GridBase>(grids.front());

    auto params_depuis_polygones = ParametresVDBDepuisMaillage();
    params_depuis_polygones.adaptrice = params->maillage_reference;
    params_depuis_polygones.bande_en_unite_globale = false;
    params_depuis_polygones.compte_voxel_bande_exterieure = static_cast<int>(largeur_de_bande);
    params_depuis_polygones.compte_voxel_bande_interieure = static_cast<int>(largeur_de_bande);
    params_depuis_polygones.genere_champs_de_distance = true;
    params_depuis_polygones.genere_volume_dense = false;
    params_depuis_polygones.utilise_grille_reference = true;
    params_depuis_polygones.grille_reference = &grille_reference_maillage;

    vdb_depuis_maillage(ctx, ctx_eval, &params_depuis_polygones, nullptr, boss, &extraction);

    if (boss.wasInterrupted()) {
        return;
    }

    auto refGrid = gridConstPtrCast<GridType>(extraction.champs_de_distance);

    auto indexGrid = extraction.grille_index;

    using BoolTreeType = typename TreeType::template ValueConverter<bool>::Type;
    typename BoolTreeType::Ptr maskTree;
    if (params->affiner_les_traits) {
        maskTree = typename BoolTreeType::Ptr(new BoolTreeType(false));
        maskTree->topologyUnion(indexGrid->tree());
        tree::LeafManager<BoolTreeType> maskLeafs(*maskTree);

        GenAdaptivityMaskOp<typename IntGridT::TreeType, BoolTreeType> op(
            AdaptriceMaillageVDB::enveloppe(params->maillage_reference),
            indexGrid->tree(),
            maskLeafs,
            params->tolerance_de_bord);
        op.run();

        tools::pruneInactive(*maskTree);

        tools::dilateActiveValues(*maskTree, 2, tools::NN_FACE, tools::IGNORE_TILES);

        mesher.setAdaptivityMask(maskTree);
    }

    if (boss.wasInterrupted()) {
        return;
    }

    mesher.setRefGrid(refGrid, params->adaptivite_interne);

    std::vector<std::string> badTransformList, badBackgroundList, badTypeList;
    for (auto &grille : grids) {
        if (boss.wasInterrupted()) {
            break;
        }

        typename GridType::ConstPtr grid = gridConstPtrCast<GridType>(grille);

        if (!grid) {
            badTypeList.push_back(grid->getName());
            continue;
        }

        if (grid->transform() != *transform) {
            badTransformList.push_back(grid->getName());
            continue;
        }

        if (!math::isApproxEqual(grid->background(), backgroundValue)) {
            badBackgroundList.push_back(grid->getName());
            continue;
        }

        mesher(*grid);

        copyMesh(maillage, mesher, grid->getName().c_str());
    }

    /* Affinage des traits. */
    if (!boss.wasInterrupted() && params->affiner_les_traits) {
        auto refGeo = AdaptriceMaillageVDB::enveloppe(params->maillage_reference);
        void *surfaceGroup = maillage.creeUnGroupeDePolygones("polygones_sur_surface");
        tbb::parallel_for(
            maillage.plagePoint(),
            SharpenFeaturesOp(
                maillage, refGeo, edgeData, *transform, surfaceGroup, maskTree.get()));
    }

    // Transfer primitive attributes
    //    if (!boss.wasInterrupted() && transferAttributes && refGeo && indexGrid) {
    //        hvdb::transferPrimitiveAttributes(*refGeo, *gdp, *indexGrid, boss, surfaceGroup);
    //    }

    if (!badTransformList.empty()) {
        std::string s = "Les grilles suivantes furent ignorées : '" +
                        outils::enchaine(badTransformList, ", ") +
                        "' car leurs transformations sont différentes de la première grille.";
        ctx_eval.rapporteAvertissement(s);
    }

    if (!badBackgroundList.empty()) {
        std::string s =
            "Les grilles suivantes furent ignorées : '" +
            outils::enchaine(badBackgroundList, ", ") +
            "' car leurs valeurs d'arrière plan sont différentes de la première grille.";
        ctx_eval.rapporteAvertissement(s);
    }

    if (!badTypeList.empty()) {
        std::string s = "Les grilles suivantes furent ignorées : '" +
                        outils::enchaine(badTypeList, ", ") +
                        "' car leurs types de données sont différents de la première grille.";
        ctx_eval.rapporteAvertissement(s);
    }
}

void maillage_depuis_vdb(ContexteKuri *ctx,
                         EnveloppeContexteEvaluation &ctx_eval,
                         ParametresVDBVersMaillage *params,
                         AdaptriceMaillage *maillage,
                         InterruptriceVDB &boss)
{
    boss.start("Conversion VDB vers polygones");

    auto grilles = outils::grilles_depuis_iteratrice(*params->groupe_grilles);
    if (grilles.empty()) {
        ctx_eval.rapporteAvertissement("Aucune grille à mailler");
        return;
    }

    auto maillage_ = AdaptriceMaillageVDB::enveloppe(maillage);

    tools::VolumeToMesh mesher(params->isovalue, params->adaptivite);
    ajoute_masque_surface(ctx_eval, params, mesher);
    ajoute_champs_adaptivite(ctx_eval, params, mesher);

    /* Crée une grille pour le maillage de référence. */
    if (params->maillage_reference) {
        // Collect all level set grids.
        std::list<GridBase::ConstPtr> grids;
        std::vector<std::string> nonLevelSetList, nonLinearList;
        for (auto grille : grilles) {
            if (boss.wasInterrupted()) {
                break;
            }

            auto &grid = grille->grid;

            const GridClass gridClass = grid->getGridClass();
            if (gridClass != GRID_LEVEL_SET) {
                nonLevelSetList.push_back(grid->getName());
                continue;
            }

            if (!grid->transform().isLinear()) {
                nonLinearList.push_back(grid->getName());
                continue;
            }

            // (We need a shallow copy to sync primitive & grid names).
            grids.push_back(grid->copyGrid());
            ConstPtrCast<GridBase>(grids.back())->setName(grid->getName());
        }

        if (!nonLevelSetList.empty()) {
            std::string s =
                "Le maillage de référence n'est supporté que pour les champs de distance "
                ", les grille suivantes furent ignorées : '" +
                outils::enchaine(nonLevelSetList, ", ") + "'.";
            ctx_eval.rapporteAvertissement(s);
        }

        if (!nonLinearList.empty()) {
            std::string s = "Les grilles suivantes furent ignorées : '" +
                            outils::enchaine(nonLinearList, ", ") +
                            "' car leurs transformations ne sont ni linéaires ni affines.";
            ctx_eval.rapporteAvertissement(s);
        }

        // Mesh using a reference surface
        if (!grids.empty() && !boss.wasInterrupted()) {
            if (grids.front()->isType<FloatGrid>()) {
                vdb_vers_polygones_reference<FloatGrid>(
                    ctx, ctx_eval, params, maillage_, mesher, grids, boss);
            }
#if 0
            // À FAIRE : double
            else if (grids.front()->isType<DoubleGrid>()) {
                vdb_vers_polygones_reference<DoubleGrid>(
                    ctx, ctx_eval, params, flux_sortie_maillage_, mesher, grids, boss);
            }
#endif
            else {
                ctx_eval.rapporteErreur("Type de grille non supporté");
            }
        }
    }
    else {
        for (auto grille : grilles) {
            if (boss.wasInterrupted()) {
                break;
            }

            if (!grille->grid) {
                continue;
            }

            grille->grid->apply<ScalarGridTypes>(mesher);
            copyMesh(maillage_, mesher, grille->grid->getName().c_str());
        }
    }

    if (boss.wasInterrupted()) {
        ctx_eval.rapporteAvertissement("Le processus fut interrompu");
    }

    boss.end();
}

}  // namespace kvdb
