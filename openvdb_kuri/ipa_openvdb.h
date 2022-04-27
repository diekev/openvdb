/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2022 Kévin Dietrich.
 * All rights reserved.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* Structure de rappels pour gérer les longs calculs. Ceci sert à interrompre au besoin lesdits
 * longs calculs. */
struct Interruptrice {
    void (*commence)(void *, const char *message);
    void (*termine)(void *);
    bool (*doit_interrompre)(void *, int pourcentage);
    void *donnees;
};

enum TypeVolume {
    INVALIDE = 0,
    R32 = 1,
    R64 = 2,
    Z32 = 3,
    Z64 = 4,
    BOOL = 5,
    VEC3_R32 = 6,
    VEC3_R64 = 7,
    VEC3_Z32 = 8,
    INDEX_POINT = 9,
    DONNEES_POINT = 10,
};

/* Structure servant à rafiner les polygones n'étant ni des triangles, ni des quadrilatères,
 * les algorithmes d'OpenVDB ne prennant pas d'autres polygones en entrée. */
struct RafineusePolygone {
    void (*ajoute_triangle)(RafineusePolygone *, long v1, long v2, long v3);
    void (*ajoute_quadrilatere)(RafineusePolygone *, long v1, long v2, long v3, long v4);

    void *donnees = nullptr;
};

struct AdaptriceMaillage {
    long (*nombre_de_polygones)(void *);
    long (*nombre_de_points)(void *);
    long (*nombre_de_sommets_polygone)(void *, long n);

    void (*point_pour_index)(void *, long n, float *x, float *y, float *z);

    void (*point_pour_sommet_polygone)(void *, long p, long s, float *x, float *y, float *z);
    void (*index_points_sommets_polygone)(void *, long n, int *index);

    /* Appelée si un polygone possède plus que 4 sommet afin que l'application cliente définissent
     * comment rafiner ces polygones. */
    void (*rafine_polygone)(void *, long n, const RafineusePolygone *);

    void *donnees = nullptr;
};

struct AdaptricePoints {
    long (*nombre_de_points)(void *);

    void (*position_pour_index)(void *, long n, float *x, float *y, float *z);

    void (*position_et_rayon_pour_index)(
        void *, long n, float *x, float *y, float *z, float *rayon);

    void (*position_rayon_et_velocite_pour_index)(void *,
                                                  long n,
                                                  float *x,
                                                  float *y,
                                                  float *z,
                                                  float *rayon,
                                                  float *vx,
                                                  float *vy,
                                                  float *vz);

    bool (*possede_attribut_pour_rayon)(void *);

    bool (*possede_attribut_pour_velocite)(void *);

    void *donnees = nullptr;
};

#ifdef __cplusplus
}
#endif
