/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#pragma once

#include "openvdb/openvdb.h"
#include "openvdb/tools/MeshToVolume.h"

struct ContexteKuri;
struct EnveloppeContexteEvaluation;
struct ExportriceGrilles;
struct InterruptriceVDB;
struct ParametresVDBDepuisMaillage;

namespace kvdb {

/* Utilisée pour que des fonctions qui appelent vdb_depuis_maillage puissent extraire
 * certaines données locales. */
struct ExtractionDonneesVDBDepuisPolygones {
    openvdb::GridBase::Ptr champs_de_distance = nullptr;
    openvdb::Int32Grid::Ptr grille_index = nullptr;
    openvdb::tools::MeshToVoxelEdgeData *donnees_aretes_voxel = nullptr;
};

void vdb_depuis_maillage(ContexteKuri *ctx,
                         EnveloppeContexteEvaluation &ctx_eval_,
                         ParametresVDBDepuisMaillage *params,
                         ExportriceGrilles *exportrice,
                         InterruptriceVDB &boss,
                         ExtractionDonneesVDBDepuisPolygones *extraction = nullptr);

}  // namespace kvdb
