/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#pragma once

struct AdaptriceSortieSpheres;
struct EnveloppeContexteEvaluation;
struct InterruptriceVDB;
struct ParametresSphereDepuisVDB;

namespace kvdb {

void sphères_depuis_vdb(EnveloppeContexteEvaluation &ctx_eval,
                        ParametresSphereDepuisVDB *params,
                        AdaptriceSortieSpheres *sortie_spheres,
                        InterruptriceVDB &boss);
}
