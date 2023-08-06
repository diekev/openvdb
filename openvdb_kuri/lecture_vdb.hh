/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#pragma once

struct ContexteKuri;
struct EnveloppeContexteEvaluation;
struct ExportriceGrilles;
struct InterruptriceVDB;
struct ParametresEcritureVDB;
struct ParametresLectureVDB;

namespace kvdb {

void lecture_vdb(ContexteKuri *ctx,
                 EnveloppeContexteEvaluation &ctx_eval,
                 ParametresLectureVDB *params,
                 ExportriceGrilles *flux_sortie_grille,
                 InterruptriceVDB &boss);

void écriture_vdb(EnveloppeContexteEvaluation &ctx_eval,
                  ParametresEcritureVDB *params,
                  InterruptriceVDB &boss);

}  // namespace kvdb
