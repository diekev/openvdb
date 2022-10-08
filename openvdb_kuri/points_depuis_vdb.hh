/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2022 Kévin Dietrich. */

#pragma once

struct ContexteKuri;
struct EnveloppeContexteEvaluation;
struct ParametresPointsDepuisVDB;
struct InterruptriceVDB;

namespace kvdb {

void points_depuis_vdb(ContexteKuri &ctx_kuri,
                       EnveloppeContexteEvaluation &ctx_eval,
                       ParametresPointsDepuisVDB &params,
                       InterruptriceVDB &boss);
}
