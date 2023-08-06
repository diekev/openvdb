/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#pragma once

struct AdaptriceMaillage;
struct ContexteKuri;
struct EnveloppeContexteEvaluation;
struct InterruptriceVDB;
struct ParametresVDBVersMaillage;

namespace kvdb {

void maillage_depuis_vdb(ContexteKuri *ctx,
                         EnveloppeContexteEvaluation &ctx_eval,
                         ParametresVDBVersMaillage *params,
                         AdaptriceMaillage *maillage,
                         InterruptriceVDB &boss);

}
