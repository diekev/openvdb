/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#pragma once

struct ContexteKuri;
struct EnveloppeContexteEvaluation;
struct AdaptricePoints;
struct ParamsVDBDepuisPoints;
struct InterruptriceVDB;

namespace kvdb {

void vdb_depuis_points(ContexteKuri *ctx,
                       EnveloppeContexteEvaluation &ctx_eval,
                       AdaptricePoints *adaptrice,
                       ParamsVDBDepuisPoints *params,
                       InterruptriceVDB &interruptrice);
}
