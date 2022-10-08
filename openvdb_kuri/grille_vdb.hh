/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2022 Kévin Dietrich. */

#pragma once

#include "openvdb/Grid.h"

struct GrilleVDB {
    openvdb::GridBase::Ptr grid;
};
