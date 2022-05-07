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

#include "ipa_openvdb.h"

#include <list>

#include "openvdb/math/Ray.h"
#include "openvdb/tools/LevelSetUtil.h"
#include "openvdb/tools/Mask.h"  // pour tools::interiorMask()
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/Morphology.h"
#include "openvdb/tools/ParticlesToLevelSet.h"
#include "openvdb/tools/VolumeToMesh.h"

#include "../InterfaceCKuri/contexte_kuri.hh"
#include "Géométrie3D/ipa.h"

#include "outils.hh"

using namespace openvdb::OPENVDB_VERSION_NAME;

/* *********************************************************************** */

// Grid type lists, for use with GEO_PrimVDB::apply(), GEOvdbApply(),
// or openvdb::GridBase::apply()

using ScalarGridTypes = openvdb::TypeList<openvdb::BoolGrid,
                                          openvdb::FloatGrid,
                                          openvdb::DoubleGrid,
                                          openvdb::Int32Grid,
                                          openvdb::Int64Grid>;

using PointGridTypes = openvdb::TypeList<openvdb::points::PointDataGrid>;

using VolumeGridTypes = ScalarGridTypes::Append<Vec3GridTypes>;

using AllGridTypes = VolumeGridTypes::Append<PointGridTypes>;

/* *********************************************************************** */

struct GrilleVDB {
    GridBase::Ptr grid;
};

GrilleVDB *VDB_copie_grille(ContexteKuri *ctx, GrilleVDB *grille)
{
    if (!grille) {
        return nullptr;
    }

    GrilleVDB *resultat = kuri_loge<GrilleVDB>(ctx);
    resultat->grid = grille->grid;
    return resultat;
}

void VDB_detruit_grille(struct ContexteKuri *ctx, struct GrilleVDB *grille)
{
    kuri_deloge(ctx, grille);
}

void VDB_accede_nom_grille(struct GrilleVDB *grille, const char **nom, long *taille)
{
    if (!grille || !grille->grid) {
        *nom = nullptr;
        *taille = 0;
        return;
    }

    const auto &nom_grille = grille->grid->getName();
    *nom = nom_grille.c_str();
    *taille = static_cast<long>(nom_grille.size());
}

void VDB_mute_nom_grille(struct GrilleVDB *grille, const char *nom, long taille)
{
    if (!grille || !grille->grid) {
        return;
    }

    const auto nouveau_nom = std::string(nom, static_cast<size_t>(taille));
    grille->grid->setName(nouveau_nom);
}

TypeVolume VDB_type_volume_pour_grille(GrilleVDB *grille)
{
    if (!grille) {
        return VOLUME_INVALIDE;
    }

    GridBase::ConstPtr grille_vdb = grille->grid;
    if (!grille_vdb) {
        return VOLUME_INVALIDE;
    }

#define VERIFIE_TYPE_ET_RETOURNE(GridType, type_volume)                                           \
    if (grille_vdb->isType<GridType>()) {                                                         \
        return type_volume;                                                                       \
    }

    VERIFIE_TYPE_ET_RETOURNE(FloatGrid, VOLUME_R32)
    VERIFIE_TYPE_ET_RETOURNE(DoubleGrid, VOLUME_R64)
    VERIFIE_TYPE_ET_RETOURNE(Int32Grid, VOLUME_Z32)
    VERIFIE_TYPE_ET_RETOURNE(Int64Grid, VOLUME_Z64)
    VERIFIE_TYPE_ET_RETOURNE(BoolGrid, VOLUME_BOOL)
    VERIFIE_TYPE_ET_RETOURNE(Vec3SGrid, VOLUME_VEC3_R32)
    VERIFIE_TYPE_ET_RETOURNE(Vec3DGrid, VOLUME_VEC3_R64)
    VERIFIE_TYPE_ET_RETOURNE(Vec3IGrid, VOLUME_VEC3_Z32)
    VERIFIE_TYPE_ET_RETOURNE(tools::PointIndexGrid, VOLUME_INDEX_POINT)
    VERIFIE_TYPE_ET_RETOURNE(points::PointDataGrid, VOLUME_DONNEES_POINT)

#undef VERIFIE_TYPE_ET_RETOURNE

    return VOLUME_INVALIDE;
}

/* *********************************************************************** */

static void exporte_grille_vdb(ContexteKuri *ctx,
                               ExportriceGrilles *exportrice,
                               GridBase::Ptr grille_vdb,
                               std::string nom)
{
    grille_vdb->setName(nom);
    GrilleVDB *grille = kuri_loge<GrilleVDB>(ctx);
    grille->grid = grille_vdb;
    exportrice->ajoute_grille(exportrice->donnees, grille);
}

#if 0
enum DrapeauxVDBDepuisMaillage {
    CHAMPS_DISTANCE_ABSOLUE = 1,
    DESACTIVE_SUPPRESSION_VOXELS_INTERSECTANT = 2,
    DESACTIVE_RENORMALISATION = 4,
    DESACTIVE_ELAGAGE_LIMITES_BANDE = 8,
};
#endif

static float calcule_taille_voxel_effective(const ParametresVDBDepuisMaillage &params)
{
    auto taille_voxel_depuis_compte = [&](float taille_axe) {
        const float dim = static_cast<float>(params.compte_voxel);
        if (params.bande_en_unite_globale) {
            const float w = params.bande_exterieure;
            return (taille_axe + 2.0f * w) / dim;
        }

        const float w = static_cast<float>(params.compte_voxel_bande_interieure);
        return taille_axe / std::max(1.0f, dim - 2.0f * w);
    };

    switch (params.methode_calcule_taille_voxel) {
        case MethodeCalculTailleVoxel::UNITE_GLOBALE:
        {
            return params.taille_voxel;
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_X:
        {
            const auto bbox = AdaptriceMaillageVDB::enveloppe(params.adaptrice).getBoundBox();
            const float taille_x = bbox.extents().x();
            return taille_voxel_depuis_compte(taille_x);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_Y:
        {
            const auto bbox = AdaptriceMaillageVDB::enveloppe(params.adaptrice).getBoundBox();
            const float taille_y = bbox.extents().y();
            return taille_voxel_depuis_compte(taille_y);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_Z:
        {
            const auto bbox = AdaptriceMaillageVDB::enveloppe(params.adaptrice).getBoundBox();
            const float taille_z = bbox.extents().z();
            return taille_voxel_depuis_compte(taille_z);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_GRAND_AXE:
        {
            const auto bbox = AdaptriceMaillageVDB::enveloppe(params.adaptrice).getBoundBox();
            const float taille_max = std::max(bbox.extents().x(),
                                              std::max(bbox.extents().y(), bbox.extents().z()));
            return taille_voxel_depuis_compte(taille_max);
        }
    }

    return 0.0f;
}

static float calcule_taille_bande_effective(const ParametresVDBDepuisMaillage &params,
                                            float bande_globale,
                                            float bande_voxel,
                                            float taille_voxel)
{
    if (params.bande_en_unite_globale) {
        return bande_globale / taille_voxel;
    }
    return bande_voxel;
}

struct DonneesCreationGrilleAttribut {
    std::string nom_attribut;
    void *poignee_attribut;
    TypeVolume type_volume;
};

static std::vector<DonneesCreationGrilleAttribut> parse_attributs_a_transferer(
    AccesseuseAttribut *accesseuse)
{
    std::vector<DonneesCreationGrilleAttribut> resultat;

    const long nombre_attributs = accesseuse->nombre_attributs(accesseuse->donnees_utilisateur);
    if (nombre_attributs == 0) {
        return resultat;
    }

    for (int i = 0; i < nombre_attributs; i++) {
        void *poignee = accesseuse->poignee_attribut(accesseuse->donnees_utilisateur, i);
        if (!poignee) {
            continue;
        }

        TypeVolume type_volume = accesseuse->type_volume_pour_attribut(poignee);
        if (type_volume == VOLUME_INVALIDE) {
            continue;
        }

        char *nom_attribut = nullptr;
        long taille_nom_attribut = 0;
        accesseuse->nom_attribut(poignee, &nom_attribut, &taille_nom_attribut);
    }

    return resultat;
}

static std::vector<Vec3s> extrait_points(const AdaptriceMaillageVDB &adaptrice_vdb,
                                         const math::Transform &transform)
{
    std::vector<Vec3s> resultat;
    const long nombre_de_points = adaptrice_vdb.pointCount();

    if (nombre_de_points == 0) {
        return resultat;
    }

    resultat.reserve(nombre_de_points);

    for (int i = 0; i < nombre_de_points; i++) {
        Vec3d point = adaptrice_vdb.getPoint(i);
        point = transform.worldToIndex(point);
        resultat.push_back(point);
    }

    return resultat;
}

struct RafineusePolygoneVDB : public RafineusePolygone {
    std::vector<Vec4I> *polygones = nullptr;
};

static std::vector<Vec4I> extrait_quads_et_triangles(const AdaptriceMaillageVDB &adaptrice_vdb)
{
    std::vector<Vec4I> resultat;
    const long nombre_de_primitives = adaptrice_vdb.polygonCount();

    if (nombre_de_primitives == 0) {
        return resultat;
    }

    resultat.reserve(nombre_de_primitives);

    RafineusePolygoneVDB rafineuse;
    rafineuse.polygones = &resultat;
    rafineuse.ajoute_triangle = [](RafineusePolygone *raf, long v0, long v1, long v2) {
        Vec4I triangle;
        triangle.x() = static_cast<int>(v0);
        triangle.y() = static_cast<int>(v1);
        triangle.z() = static_cast<int>(v2);
        triangle.w() = util::INVALID_IDX;
        static_cast<RafineusePolygoneVDB *>(raf)->polygones->push_back(triangle);
    };
    rafineuse.ajoute_quadrilatere =
        [](RafineusePolygone *raf, long v0, long v1, long v2, long v3) {
            Vec4I triangle;
            triangle.x() = static_cast<int>(v0);
            triangle.y() = static_cast<int>(v1);
            triangle.z() = static_cast<int>(v2);
            triangle.w() = static_cast<int>(v3);
            static_cast<RafineusePolygoneVDB *>(raf)->polygones->push_back(triangle);
        };

    for (int i = 0; i < nombre_de_primitives; i++) {
        const long nombre_de_sommets = adaptrice_vdb.vertexCount(i);

        if (nombre_de_sommets == 3) {
            int index[3];
            adaptrice_vdb.indexSommetsPolygone(i, index);
            Vec4I triangle;
            triangle.x() = index[0];
            triangle.y() = index[1];
            triangle.z() = index[2];
            triangle.w() = util::INVALID_IDX;
            resultat.push_back(triangle);
        }
        else if (nombre_de_sommets == 4) {
            int index[4];
            adaptrice_vdb.indexSommetsPolygone(i, index);
            Vec4I quad;
            quad.x() = index[0];
            quad.y() = index[1];
            quad.z() = index[2];
            quad.w() = index[3];
            resultat.push_back(quad);
        }
        else if (nombre_de_sommets > 4) {
            adaptrice_vdb.rafinePolygone(i, rafineuse);
        }
    }

    return resultat;
}

/* Utilisée pour que des fonctions qui appelent VDB_depuis_polygones_impl puissent extraire
 * certaines données locales. */
struct ExtractionDonneesVDBDepuisPolygones {
    GridBase::Ptr champs_de_distance = nullptr;
    Int32Grid::Ptr grille_index = nullptr;
    tools::MeshToVoxelEdgeData *donnees_aretes_voxel = nullptr;
};

static void VDB_depuis_polygones_impl(ContexteKuri *ctx,
                                      EnveloppeContexteEvaluation &ctx_eval_,
                                      ParametresVDBDepuisMaillage *params,
                                      ExportriceGrilles *exportrice,
                                      InterruptriceVDB &boss,
                                      ExtractionDonneesVDBDepuisPolygones *extraction = nullptr)
{
    boss.start("Conversion de polygones en volumes VDB");

    if (!params->genere_champs_de_distance && !params->genere_volume_dense &&
        !params->transfere_attributs) {
        ctx_eval_.rapporteErreur("Rien à générer");
        return;
    }

    math::Transform::Ptr transform;
    float taille_voxel;

    if (params->utilise_grille_reference) {
        if (!params->grille_reference || !params->grille_reference->grid) {
            ctx_eval_.rapporteErreur("Aucune grille de reference trouvée");
            return;
        }

        transform = params->grille_reference->grid->transform().copy();
        taille_voxel = static_cast<float>(transform->voxelSize()[0]);
    }
    else {
        taille_voxel = calcule_taille_voxel_effective(*params);
        transform = math::Transform::createLinearTransform(taille_voxel);
    }

    if (taille_voxel < 1e-5) {
        ctx_eval_.rapporteErreur("La taille de voxel (", taille_voxel, ") est trop petite");
        return;
    }

    AdaptriceMaillageVDB adaptrice_vdb = AdaptriceMaillageVDB::enveloppe(params->adaptrice);

    std::vector<Vec3s> points = extrait_points(adaptrice_vdb, *transform);
    std::vector<Vec4I> polygones = extrait_quads_et_triangles(adaptrice_vdb);

    /* Même si nous n'avons aucun point ou polygone, nous créons les grilles, car des opérations
     * suivantes pourraient tout de même gérer des grilles vides. */
    if (points.empty()) {
        ctx_eval_.rapporteAvertissement("Aucun point dans le maillage d'entrée");
    }
    if (polygones.empty()) {
        ctx_eval_.rapporteAvertissement("Aucun polygone à convertir");
    }

    tools::QuadAndTriangleDataAdapter<Vec3s, Vec4I> mesh(points, polygones);

    /* Grille d'index pour les attributs. */
    Int32Grid::Ptr grille_index;
    if (params->transfere_attributs || extraction) {
        grille_index.reset(new Int32Grid(0));
    }

    if (extraction) {
        extraction->grille_index = grille_index;
        if (extraction->donnees_aretes_voxel) {
            extraction->donnees_aretes_voxel->convert(points, polygones);
        }
    }

    int drapeaux = params->champs_distance_absolue ? tools::UNSIGNED_DISTANCE_FIELD : 0;

    const float bande_exterieure = calcule_taille_bande_effective(
        *params, params->bande_exterieure, params->compte_voxel_bande_exterieure, taille_voxel);
    const float bande_interieure = params->remplis_interieur ?
                                       std::numeric_limits<float>::max() :
                                       calcule_taille_bande_effective(
                                           *params,
                                           params->bande_interieure,
                                           params->compte_voxel_bande_interieure,
                                           taille_voxel);

    FloatGrid::Ptr grille = tools::meshToVolume<FloatGrid>(
        boss, mesh, *transform, bande_exterieure, bande_interieure, drapeaux, grille_index.get());

    if (extraction) {
        extraction->champs_de_distance = grille;
    }

    if (!boss.wasInterrupted() && params->genere_champs_de_distance && exportrice) {
        exporte_grille_vdb(ctx,
                           exportrice,
                           grille,
                           outils::chaine_depuis_accesseuse(params->nom_champs_distance));
    }

    if (params->genere_volume_dense && params->champs_distance_absolue) {
        ctx_eval_.rapporteAvertissement("Impossible de générer un volume dense car le champs de "
                                        "distance à générer possède des valeurs absolues");
    }

    if (!boss.wasInterrupted() && params->genere_volume_dense &&
        !params->champs_distance_absolue) {
        /* Si nous n'exportons pas le champs de distance, nous pouvons le modifier directemet. */
        FloatGrid::Ptr grille_fog;
        if (params->genere_champs_de_distance) {
            grille_fog = grille->deepCopy();
        }
        else {
            grille_fog = grille;
        }

        tools::sdfToFogVolume(*grille_fog);

        exporte_grille_vdb(ctx,
                           exportrice,
                           grille_fog,
                           outils::chaine_depuis_accesseuse(params->nom_volume_dense));
    }

    if (!boss.wasInterrupted() && params->transfere_attributs) {
        // À FAIRE
    }

    if (boss.wasInterrupted()) {
        ctx_eval_.rapporteAvertissement("Le processus fut interrompu");
    }

    boss.end();
}

void VDB_depuis_polygones(ContexteKuri *ctx,
                          ContexteEvaluation *ctx_eval,
                          ParametresVDBDepuisMaillage *params,
                          ExportriceGrilles *exportrice,
                          Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        VDB_depuis_polygones_impl(ctx, ctx_eval_, params, exportrice, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

static std::vector<GrilleVDB *> grilles_depuis_iteratrice(IteratriceGrillesVDB &iteratrice)
{
    std::vector<GrilleVDB *> resultat;

    auto grille = iteratrice.suivante(iteratrice.donnees_utilisateur);
    while (grille) {
        resultat.push_back(grille);
        grille = iteratrice.suivante(iteratrice.donnees_utilisateur);
    }

    return resultat;
}

struct InteriorMaskOp {
    InteriorMaskOp(double iso = 0.0) : inIsovalue(iso)
    {
    }

    template <typename GridType>
    void operator()(const GridType &grid)
    {
        outGridPtr = openvdb::tools::interiorMask(grid, inIsovalue);
    }

    const double inIsovalue;
    openvdb::BoolGrid::Ptr outGridPtr;
};

// Extract a boolean mask from a grid of any type.
inline GridBase::ConstPtr getMaskFromGrid(const GridBase::ConstPtr &gridPtr, double isovalue = 0.0)
{
    if (!gridPtr) {
        return nullptr;
    }

    if (gridPtr->isType<openvdb::BoolGrid>()) {
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
        ctx_eval.rapporteAvertissement("Aucune grille pour le masque de surface");
        return;
    }

    if (VDB_type_volume_pour_grille(params->grille_champs_adaptivite) != VOLUME_R32) {
        ctx_eval.rapporteAvertissement("La grille du champs d'adaptivité n'est pas de type réel");
        return;
    }

    auto grille = gridConstPtrCast<FloatGrid>(params->grille_champs_adaptivite->grid);
    mesher.setSpatialAdaptivity(grille);
}

void copyMesh(AdaptriceMaillageVDB &maillage, tools::VolumeToMesh &mesher, const char *gridName)
{
    /* Exporte les points. */
    const openvdb::tools::PointList &points = mesher.pointList();

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
    openvdb::tools::PolygonPoolList &polygonPoolList = mesher.polygonPoolList();

    const char exteriorFlag = char(openvdb::tools::POLYFLAG_EXTERIOR);
    const char seamLineFlag = char(openvdb::tools::POLYFLAG_FRACTURE_SEAM);

    // index 0 --> interior, not on seam
    // index 1 --> interior, on seam
    // index 2 --> surface,  not on seam
    // index 3 --> surface,  on seam
    long nquads[4] = {0, 0, 0, 0};
    long ntris[4] = {0, 0, 0, 0};
    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
        const openvdb::tools::PolygonPool &polygons = polygonPoolList[n];
        for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
            int flags = (((polygons.quadFlags(i) & exteriorFlag) != 0) << 1) |
                        ((polygons.quadFlags(i) & seamLineFlag) != 0);
            ++nquads[flags];
        }
        for (size_t i = 0, I = polygons.numTriangles(); i < I; ++i) {
            int flags = (((polygons.triangleFlags(i) & exteriorFlag) != 0) << 1) |
                        ((polygons.triangleFlags(i) & seamLineFlag) != 0);
            ++ntris[flags];
        }
    }

    long nverts[4] = {nquads[0] * 4 + ntris[0] * 3,
                      nquads[1] * 4 + ntris[1] * 3,
                      nquads[2] * 4 + ntris[2] * 3,
                      nquads[3] * 4 + ntris[3] * 3};
    std::vector<int> verts[4];
    for (int flags = 0; flags < 4; ++flags) {
        verts[flags].resize(nverts[flags]);
    }

    long iquad[4] = {0, 0, 0, 0};
    long itri[4] = {nquads[0] * 4, nquads[1] * 4, nquads[2] * 4, nquads[3] * 4};

    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
        const openvdb::tools::PolygonPool &polygons = polygonPoolList[n];

        // Copy quads
        for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
            const openvdb::Vec4I &quad = polygons.quad(i);
            int flags = (((polygons.quadFlags(i) & exteriorFlag) != 0) << 1) |
                        ((polygons.quadFlags(i) & seamLineFlag) != 0);
            verts[flags][iquad[flags]++] = quad[0] + decalage_point;
            verts[flags][iquad[flags]++] = quad[1] + decalage_point;
            verts[flags][iquad[flags]++] = quad[2] + decalage_point;
            verts[flags][iquad[flags]++] = quad[3] + decalage_point;
        }

        // Copy triangles (adaptive mesh)
        for (size_t i = 0, I = polygons.numTriangles(); i < I; ++i) {
            const openvdb::Vec3I &triangle = polygons.triangle(i);
            int flags = (((polygons.triangleFlags(i) & exteriorFlag) != 0) << 1) |
                        ((polygons.triangleFlags(i) & seamLineFlag) != 0);
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
    using BoolLeafManager = openvdb::tree::LeafManager<BoolTreeType>;

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
    using IndexAccessorType = typename openvdb::tree::ValueAccessor<const IndexTreeType>;
    IndexAccessorType idxAcc(mIndexTree);

    Vec3s tmpN, normal;
    int tmpIdx;

    openvdb::Coord ijk, nijk;
    typename BoolTreeType::LeafNodeType::ValueOnIter iter;

    for (size_t n = range.begin(); n < range.end(); ++n) {
        iter = mLeafs.leaf(n).beginValueOn();
        for (; iter; ++iter) {
            ijk = iter.getCoord();

            bool edgeVoxel = false;

            int idx = idxAcc.getValue(ijk);

            normal = normaux[idx];

            for (size_t i = 0; i < 18; ++i) {
                nijk = ijk + openvdb::util::COORD_OFFSETS[i];
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
    using EdgeData = openvdb::tools::MeshToVoxelEdgeData;

    SharpenFeaturesOp(AdaptriceMaillageVDB &meshGeo,
                      const AdaptriceMaillageVDB &refGeo,
                      EdgeData &edgeData,
                      const openvdb::math::Transform &xform,
                      const void *surfacePrims = nullptr,
                      const openvdb::BoolTree *mask = nullptr);

    void operator()(const tbb::blocked_range<long> &) const;

  private:
    AdaptriceMaillageVDB &mMeshGeo;
    const AdaptriceMaillageVDB &mRefGeo;
    EdgeData &mEdgeData;
    const openvdb::math::Transform &mXForm;
    const void *mSurfacePrims;
    const openvdb::BoolTree *mMaskTree;
};

SharpenFeaturesOp::SharpenFeaturesOp(AdaptriceMaillageVDB &meshGeo,
                                     const AdaptriceMaillageVDB &refGeo,
                                     EdgeData &edgeData,
                                     const openvdb::math::Transform &xform,
                                     const void *surfacePrims,
                                     const openvdb::BoolTree *mask)
    : mMeshGeo(meshGeo), mRefGeo(refGeo), mEdgeData(edgeData), mXForm(xform),
      mSurfacePrims(surfacePrims), mMaskTree(mask)
{
}

void SharpenFeaturesOp::operator()(const tbb::blocked_range<long> &range) const
{
    openvdb::tools::MeshToVoxelEdgeData::Accessor acc = mEdgeData.getAccessor();

    using BoolAccessor = openvdb::tree::ValueAccessor<const openvdb::BoolTree>;
    std::unique_ptr<BoolAccessor> maskAcc;

    if (mMaskTree) {
        maskAcc.reset(new BoolAccessor(*mMaskTree));
    }

    Vec3s tmpN, tmpP, avgP;
    math::BBox<Vec3d> cell;

    openvdb::Vec3d pos, normal;
    openvdb::Coord ijk;

    std::vector<openvdb::Vec3d> points, normals;
    std::vector<openvdb::Index32> primitives;

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

        pos = openvdb::tools::findFeaturePoint(points, normals);

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
                                         std::list<openvdb::GridBase::ConstPtr> const &grids,
                                         InterruptriceVDB &boss)
{
    typename GridType::ConstPtr firstGrid = openvdb::gridConstPtrCast<GridType>(grids.front());

    if (!firstGrid) {
        ctx_eval.rapporteErreur("Type de grille non supporté");
        return;
    }

    using TreeType = typename GridType::TreeType;
    using ValueType = typename GridType::ValueType;
    using IntGridT = typename GridType::template ValueConverter<openvdb::Int32>::Type;
    openvdb::math::Transform::Ptr transform = firstGrid->transform().copy();
    const ValueType backgroundValue = firstGrid->background();
    const openvdb::GridClass gridClass = firstGrid->getGridClass();

    /* Crée la grille de référence. */
    openvdb::tools::MeshToVoxelEdgeData edgeData;

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

    VDB_depuis_polygones_impl(ctx, ctx_eval, &params_depuis_polygones, nullptr, boss, &extraction);

    if (boss.wasInterrupted()) {
        return;
    }

    auto refGrid = openvdb::gridConstPtrCast<GridType>(extraction.champs_de_distance);

    auto indexGrid = extraction.grille_index;

    using BoolTreeType = typename TreeType::template ValueConverter<bool>::Type;
    typename BoolTreeType::Ptr maskTree;
    if (params->affiner_les_traits) {
        maskTree = typename BoolTreeType::Ptr(new BoolTreeType(false));
        maskTree->topologyUnion(indexGrid->tree());
        openvdb::tree::LeafManager<BoolTreeType> maskLeafs(*maskTree);

        GenAdaptivityMaskOp<typename IntGridT::TreeType, BoolTreeType> op(
            AdaptriceMaillageVDB::enveloppe(params->maillage_reference),
            indexGrid->tree(),
            maskLeafs,
            params->tolerance_de_bord);
        op.run();

        openvdb::tools::pruneInactive(*maskTree);

        openvdb::tools::dilateActiveValues(
            *maskTree, 2, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);

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

        typename GridType::ConstPtr grid = openvdb::gridConstPtrCast<GridType>(grille);

        if (!grid) {
            badTypeList.push_back(grid->getName());
            continue;
        }

        if (grid->transform() != *transform) {
            badTransformList.push_back(grid->getName());
            continue;
        }

        if (!openvdb::math::isApproxEqual(grid->background(), backgroundValue)) {
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

void VDB_vers_polygones_impl(ContexteKuri *ctx,
                             EnveloppeContexteEvaluation &ctx_eval,
                             ParametresVDBVersMaillage *params,
                             AdaptriceMaillage *maillage,
                             InterruptriceVDB &boss)
{
    boss.start("Conversion VDB vers polygones");

    auto grilles = grilles_depuis_iteratrice(*params->groupe_grilles);
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
        std::list<openvdb::GridBase::ConstPtr> grids;
        std::vector<std::string> nonLevelSetList, nonLinearList;
        for (auto grille : grilles) {
            if (boss.wasInterrupted()) {
                break;
            }

            auto &grid = grille->grid;

            const openvdb::GridClass gridClass = grid->getGridClass();
            if (gridClass != openvdb::GRID_LEVEL_SET) {
                nonLevelSetList.push_back(grid->getName());
                continue;
            }

            if (!grid->transform().isLinear()) {
                nonLinearList.push_back(grid->getName());
                continue;
            }

            // (We need a shallow copy to sync primitive & grid names).
            grids.push_back(grid->copyGrid());
            openvdb::ConstPtrCast<openvdb::GridBase>(grids.back())->setName(grid->getName());
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
            if (grids.front()->isType<openvdb::FloatGrid>()) {
                vdb_vers_polygones_reference<openvdb::FloatGrid>(
                    ctx, ctx_eval, params, maillage_, mesher, grids, boss);
            }
#if 0
            // À FAIRE : double
            else if (grids.front()->isType<openvdb::DoubleGrid>()) {
                vdb_vers_polygones_reference<openvdb::DoubleGrid>(
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

void VDB_vers_polygones(ContexteKuri *ctx,
                        ContexteEvaluation *ctx_eval,
                        ParametresVDBVersMaillage *params,
                        AdaptriceMaillage *maillage,
                        Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        VDB_vers_polygones_impl(ctx, ctx_eval_, params, maillage, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

struct ListeParticules : public AdaptricePoints {
    // Interface de convenience.
    bool possede_rayon() const
    {
        return this->possede_attribut_pour_rayon &&
               this->possede_attribut_pour_rayon(this->donnees);
    }

    bool possede_velocite() const
    {
        return this->possede_attribut_pour_velocite &&
               this->possede_attribut_pour_velocite(this->donnees);
    }

    // Interface requise par OpenVDB.
    size_t size() const
    {
        return static_cast<size_t>(this->nombre_de_points(this->donnees));
    }

    void getPos(size_t n, Vec3R &xyz) const
    {
        float x, y, z;
        this->position_pour_index(this->donnees, static_cast<size_t>(n), &x, &y, &z);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
    }

    void getPosRad(size_t n, Vec3R &xyz, Real &radius) const
    {
        float x, y, z, rayon;
        this->position_et_rayon_pour_index(
            this->donnees, static_cast<size_t>(n), &x, &y, &z, &rayon);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
        radius = rayon;
    }

    void getPosRadVel(size_t n, Vec3R &xyz, Real &radius, Vec3R &velocity) const
    {
        float x, y, z, rayon, vx, vy, vz;
        this->position_rayon_et_velocite_pour_index(
            this->donnees, static_cast<size_t>(n), &x, &y, &z, &rayon, &vx, &vy, &vz);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
        radius = rayon;
        velocity.x() = vx;
        velocity.y() = vy;
        velocity.z() = vz;
    }

    /* Pour le transfert d'attributs, nous utilisons une grille d'index qui sera ensuite utilisée
     * pour transférer les attributs. Cette grille stocke l'index du point le plus proche de chaque
     * voxel. */
    void getAtt(size_t n, Int32 &att) const
    {
        att = static_cast<Int32>(n);
    }
};

struct ParamsVDBDepuisParticules {
    float poids_rayon;
    float poids_velocity;
    bool cree_des_trainees;
    bool transfere_attributs;

    AccesseuseAttribut *acces_attributs_points;
};

void VDB_depuis_particules(AdaptricePoints *adaptrice,
                           ParamsVDBDepuisParticules *params,
                           Interruptrice *interruptrice)
{
    //    ListeParticules liste_particules = ListeParticules{*adaptrice};
    //    InterruptriceVDB boss{interruptrice};

    //    tools::ParticlesToLevelSet p;
}

static void VDB_depuis_fichier_impl(ContexteKuri *ctx,
                                    EnveloppeContexteEvaluation &ctx_eval,
                                    ParametresLectureVDB *params,
                                    ExportriceGrilles *flux_sortie_grille,
                                    InterruptriceVDB &boss)
{
    boss.start("Lecture d'un fichier OpenVDB");

    const std::string chemin_fichier = outils::chaine_depuis_accesseuse(params->chemin_fichier);
    const Index64 limite_copie = static_cast<Index64>(1.0e9 * params->limite_pour_copier);

    io::File fichier(chemin_fichier);
    fichier.setCopyMaxBytes(limite_copie);
    fichier.open(params->chargement_tardif);

    MetaMap::Ptr metadonnees_fichier = fichier.getMetadata();
    if (!metadonnees_fichier) {
        metadonnees_fichier.reset(new MetaMap);
    }

    bool rogne = params->rogne;
    BBoxd boite_rognage;
    if (rogne && params->maillage_reference) {
        AdaptriceMaillageVDB maillage = AdaptriceMaillageVDB::enveloppe(
            params->maillage_reference);
        boite_rognage = maillage.getBoundBox();
        rogne = boite_rognage.isSorted();
    }

    for (io::File::NameIterator iter_nom = fichier.beginName(); iter_nom != fichier.endName();
         ++iter_nom) {
        if (boss.wasInterrupted()) {
            throw std::runtime_error("Le processus fut interrompu");
        }

        const std::string &nom_grille = iter_nom.gridName();

        GridBase::Ptr grille;
        if (params->metadonnees_seules) {
            grille = fichier.readGridMetadata(nom_grille);
        }
        else if (rogne) {
            grille = fichier.readGrid(nom_grille, boite_rognage);
        }
        else {
            grille = fichier.readGrid(nom_grille);
        }

        if (!grille) {
            continue;
        }

        /* Transfère les métadonnées du fichier à la grille. */
        for (MetaMap::ConstMetaIterator iter_metadonnee = metadonnees_fichier->beginMeta(),
                                        end = metadonnees_fichier->endMeta();
             iter_metadonnee != end;
             ++iter_metadonnee) {
            /* Ne transfère une métadonnée que si la grille n'en a pas une avec un nom similaire.
             */
            if (Metadata::Ptr meta = iter_metadonnee->second) {
                const std::string name = iter_metadonnee->first;
                if (!(*grille)[name]) {
                    grille->insertMeta(name, *meta);
                }
            }
        }

        exporte_grille_vdb(ctx, flux_sortie_grille, grille, nom_grille);
    }
    fichier.close();

    boss.end();
}

void VDB_depuis_fichier(struct ContexteKuri *ctx,
                        struct ContexteEvaluation *ctx_eval,
                        struct ParametresLectureVDB *params,
                        struct ExportriceGrilles *flux_sortie_grille,
                        struct Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        VDB_depuis_fichier_impl(ctx, ctx_eval_, params, flux_sortie_grille, boss);
    }
    catch (IoError &e) {
        std::string message = e.what();
        /* Enlève le "IOError: " du message. */
        message = message.substr(9);

        if (params->comportement_fichier_manquant ==
            ComportementFichierManquant::RAPPORTE_ERREUR) {
            ctx_eval_.rapporteErreur(message);
        }
        else {
            ctx_eval_.rapporteAvertissement(message);
        }
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}
