/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "vdb_depuis_maillage.hh"

#include "openvdb/tools/LevelSetUtil.h"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"
#include "outils.hh"

using namespace openvdb;

namespace kvdb {

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

void vdb_depuis_maillage(ContexteKuri *ctx,
                         EnveloppeContexteEvaluation &ctx_eval_,
                         ParametresVDBDepuisMaillage *params,
                         ExportriceGrilles *exportrice,
                         InterruptriceVDB &boss,
                         ExtractionDonneesVDBDepuisPolygones *extraction)
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
        outils::exporte_grille_vdb(ctx,
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

        outils::exporte_grille_vdb(ctx,
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

}  // namespace kvdb
