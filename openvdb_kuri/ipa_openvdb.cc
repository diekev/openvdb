/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2022 Kévin Dietrich. */

#include "ipa_openvdb.h"

#include <list>

#include "openvdb/math/Ray.h"
#include "openvdb/tools/Dense.h"
#include "openvdb/tools/LevelSetUtil.h"
#include "openvdb/tools/Mask.h"  // pour tools::interiorMask()
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/Morphology.h"
#include "openvdb/tools/ParticlesToLevelSet.h"
#include "openvdb/tools/VolumeToMesh.h"

#include "../InterfaceCKuri/contexte_kuri.hh"
#include "Géométrie3D/ipa.h"

#include "grille_vdb.hh"
#include "lecture_vdb.hh"
#include "maillage_depuis_vdb.hh"
#include "outils.hh"
#include "points_depuis_vdb.hh"
#include "vdb_depuis_maillage.hh"
#include "vdb_depuis_points.hh"

using namespace openvdb;

/* ------------------------------------------------------------------------- */
/** \name Utilitaires.
 * \{ */

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

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Nommage d'une grille.
 * \{ */

void VDB_donne_nom_grille(struct GrilleVDB *grille, const char **nom, long *taille)
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

void VDB_definis_nom_grille(struct GrilleVDB *grille, const char *nom, long taille)
{
    if (!grille || !grille->grid) {
        return;
    }

    const auto nouveau_nom = std::string(nom, static_cast<size_t>(taille));
    grille->grid->setName(nouveau_nom);
}

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Dimension mondiale d'une grille.
 * \{ */

void VDB_donne_dimension_monde(GrilleVDB *grille, DimensionMondeVDB *r_dimension)
{
    if (!grille || !grille->grid || !r_dimension) {
        if (r_dimension) {
            r_dimension->min_x = 0.0f;
            r_dimension->min_y = 0.0f;
            r_dimension->min_z = 0.0f;
            r_dimension->max_x = 0.0f;
            r_dimension->max_y = 0.0f;
            r_dimension->max_z = 0.0f;
        }
        return;
    }

    auto const &grid = grille->grid;
    auto const &xform = grid->transform();

    auto const bbox = grid->evalActiveVoxelBoundingBox();
    auto const min_monde = xform.indexToWorld(bbox.min());
    auto const max_monde = xform.indexToWorld(bbox.max());

    r_dimension->min_x = min_monde.x();
    r_dimension->min_y = min_monde.y();
    r_dimension->min_z = min_monde.z();

    r_dimension->max_x = max_monde.x();
    r_dimension->max_y = max_monde.y();
    r_dimension->max_z = max_monde.z();
}

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Conversion vers tampon dense.
 * \{ */

void VDB_donne_tampon_dense(ContexteKuri *ctx,
                            GrilleVDB *grille,
                            DonneesConversionVersDense *r_donnees)
{
    if (!grille || !grille->grid || !r_donnees) {
        if (r_donnees) {
            r_donnees->donnees = nullptr;
            r_donnees->taille_donnees = 0;
            r_donnees->dim_x = 0;
            r_donnees->dim_y = 0;
            r_donnees->dim_z = 0;
        }
        return;
    }

    auto const &grid = grille->grid;
    auto bbox = grid->evalActiveVoxelBoundingBox();

    float *r_data = kuri_loge_tableau<float>(ctx, bbox.volume());

    auto dense = openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ>(bbox, r_data);

    openvdb::tools::copyToDense(*openvdb::gridConstPtrCast<openvdb::FloatGrid>(grid), dense);

    r_donnees->donnees = r_data;
    r_donnees->taille_donnees = bbox.volume();
    r_donnees->dim_x = bbox.dim().x();
    r_donnees->dim_y = bbox.dim().y();
    r_donnees->dim_z = bbox.dim().z();
}

void VDB_detruit_tampon_dense(ContexteKuri *ctx, DonneesConversionVersDense *donnees)
{
    if (!donnees) {
        return;
    }

    kuri_déloge_tableau(ctx, donnees->donnees, donnees->taille_donnees);
}

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Conversions VDB <-> maillage.
 * \{ */

void VDB_depuis_polygones(ContexteKuri *ctx,
                          ContexteEvaluation *ctx_eval,
                          ParametresVDBDepuisMaillage *params,
                          ExportriceGrilles *exportrice,
                          Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        kvdb::vdb_depuis_maillage(ctx, ctx_eval_, params, exportrice, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
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
        kvdb::maillage_depuis_vdb(ctx, ctx_eval_, params, maillage, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Conversions VDB <-> points.
 * \{ */

void VDB_depuis_points(ContexteKuri *ctx,
                       ContexteEvaluation *ctx_eval,
                       AdaptricePoints *adaptrice,
                       ParamsVDBDepuisPoints *params,
                       Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        kvdb::vdb_depuis_points(ctx, ctx_eval_, adaptrice, params, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

void VDB_distribue_points(ContexteKuri *ctx,
                          ContexteEvaluation *ctx_eval,
                          ParametresPointsDepuisVDB *params,
                          Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        kvdb::points_depuis_vdb(*ctx, ctx_eval_, *params, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Lecture fichier.
 * \{ */

void VDB_depuis_fichier(struct ContexteKuri *ctx,
                        struct ContexteEvaluation *ctx_eval,
                        struct ParametresLectureVDB *params,
                        struct ExportriceGrilles *flux_sortie_grille,
                        struct Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluation ctx_eval_ = EnveloppeContexteEvaluation::enveloppe(ctx_eval);

    try {
        kvdb::lecture_vdb(ctx, ctx_eval_, params, flux_sortie_grille, boss);
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

/** \} */
