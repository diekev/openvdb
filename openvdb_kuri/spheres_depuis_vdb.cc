/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "spheres_depuis_vdb.hh"

#include "openvdb/tools/VolumeToSpheres.h"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"
#include "outils.hh"

using namespace openvdb;

namespace kvdb {

static Vec2s détermine_limites_rayon(ParametresSphereDepuisVDB *params)
{
    float min = std::numeric_limits<float>::min();
    float max = std::numeric_limits<float>::max();

    if (params->utilise_rayon_minimum) {
        min = params->rayon_minimum;
    }

    if (params->utilise_rayon_maximum) {
        max = params->rayon_maximum;
    }

    return Vec2s(min, max);
}

static Vec2i détermine_compte_de_sphères(ParametresSphereDepuisVDB *params)
{
    int min = 0;
    int max = params->compte_de_points;

    if (params->utilise_nombre_minimum) {
        min = params->nombre_minimum;
    }
    if (params->utilise_nombre_maximum) {
        max = params->nombre_maximum;
    }
    return Vec2i(min, max);
}

static void rapporte_avertissement_grilles_ignorées(
    EnveloppeContexteEvaluation &ctx_eval,
    std::vector<std::string> const &grilles_ignorées,
    std::string const &message)
{
    if (!grilles_ignorées.empty()) {
        return;
    }

    std::stringstream ss;
    ss << message;

    ss << " Les grilles suivantes furent ignorées : ";

    auto virgule = "'";

    for (auto const &nom : grilles_ignorées) {
        ss << virgule << nom;
        virgule = ", ";
    }

    ss << "'";

    ctx_eval.rapporteAvertissement(ss.str());
}

static void rapporte_interruption(EnveloppeContexteEvaluation &ctx_eval, InterruptriceVDB &boss)
{
    if (boss.wasInterrupted()) {
        ctx_eval.rapporteAvertissement("Processus interrompu");
    }
}

void sphères_depuis_vdb(EnveloppeContexteEvaluation &ctx_eval,
                        ParametresSphereDepuisVDB *params,
                        AdaptriceSortieSpheres *sortie_spheres,
                        InterruptriceVDB &boss)
{
    auto grilles = outils::grilles_depuis_iteratrice(*params->grilles);
    if (grilles.empty()) {
        ctx_eval.rapporteAvertissement("Aucune grille trouvée.");
        return;
    }

    auto const limites_rayon = détermine_limites_rayon(params);
    auto const compte_de_sphères = détermine_compte_de_sphères(params);

    std::vector<std::string> grilles_ignorées;

    int identifiant_grille = 1;
    for (auto grille : grilles) {
        if (boss.wasInterrupted()) {
            break;
        }

        Vec2s limites_rayon_grille = limites_rayon;
        if (params->utilise_unites_mondiales) {
            auto const taille_voxel = float(grille->grid->voxelSize()[0]);
            limites_rayon_grille *= taille_voxel;
        }

        limites_rayon_grille[1] = std::max(limites_rayon_grille[1],
                                           limites_rayon_grille[0] + float(1e-5));

        std::vector<openvdb::Vec4s> sphères;
        if (grille->grid->type() == FloatGrid::gridType()) {
            FloatGrid::ConstPtr gridPtr = gridConstPtrCast<FloatGrid>(grille->grid);
            tools::fillWithSpheres(*gridPtr,
                                   sphères,
                                   compte_de_sphères,
                                   params->permet_superposition,
                                   limites_rayon_grille[0],
                                   limites_rayon_grille[1],
                                   params->valeur_isometrique,
                                   params->compte_de_points,
                                   &boss);
        }
        else if (grille->grid->type() == DoubleGrid::gridType()) {
            DoubleGrid::ConstPtr gridPtr = gridConstPtrCast<DoubleGrid>(grille->grid);
            tools::fillWithSpheres(*gridPtr,
                                   sphères,
                                   compte_de_sphères,
                                   params->permet_superposition,
                                   limites_rayon_grille[0],
                                   limites_rayon_grille[1],
                                   params->valeur_isometrique,
                                   params->compte_de_points,
                                   &boss);
        }
        else {
            grilles_ignorées.push_back(grille->grid->getName());
            continue;
        }

        /* Copie les sphères. */
        for (auto sphère : sphères) {
            auto index_sphère = sortie_spheres->ajoute_sphere(
                sortie_spheres, sphère[0], sphère[1], sphère[2], sphère[3]);

            if (params->ajoute_attribut_id) {
                sortie_spheres->definis_attribut_identifiant(
                    sortie_spheres, index_sphère, identifiant_grille);
            }

            if (params->ajoute_attribut_rayon) {
                sortie_spheres->definis_attribut_rayon(sortie_spheres, index_sphère, sphère[3]);
            }
        }
    }

    rapporte_avertissement_grilles_ignorées(
        ctx_eval,
        grilles_ignorées,
        "Seules les grilles scalaires (float/double) sont supportées.");

    rapporte_interruption(ctx_eval, boss);
}

}  // namespace kvdb
