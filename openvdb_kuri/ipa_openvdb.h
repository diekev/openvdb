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
#else
typedef char bool;
#endif

struct ContexteKuri;

struct GrilleVDB;

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

void VDB_detruit_grille(struct ContexteKuri *ctx, struct GrilleVDB *grille);
void VDB_accede_nom_grille(struct GrilleVDB *grille, const char **nom, long *taille);
void VDB_mute_nom_grille(struct GrilleVDB *grille, const char *nom, long taille);
enum TypeVolume VDB_type_volume_pour_grille(struct GrilleVDB *grille);

struct ExportriceGrilles {
    void (*ajoute_grille)(void *, struct GrilleVDB *);
    void *donnees;
};

/* Structure de rappels pour gérer les longs calculs. Ceci sert à interrompre au besoin lesdits
 * longs calculs. */
struct Interruptrice {
    void (*commence)(void *, const char *message);
    void (*termine)(void *);
    bool (*doit_interrompre)(void *, int pourcentage);
    void *donnees;
};

/* Structure servant à rafiner les polygones n'étant ni des triangles, ni des quadrilatères,
 * les algorithmes d'OpenVDB ne prennant pas d'autres polygones en entrée. */
struct RafineusePolygone {
    void (*ajoute_triangle)(struct RafineusePolygone *, long v1, long v2, long v3);
    void (*ajoute_quadrilatere)(struct RafineusePolygone *, long v1, long v2, long v3, long v4);

    void *donnees;
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
    void (*rafine_polygone)(void *, long n, const struct RafineusePolygone *);

    void (*calcule_boite_englobante)(void *,
                                     float *min_x,
                                     float *min_y,
                                     float *min_z,
                                     float *max_x,
                                     float *max_y,
                                     float *max_z);

    void *donnees;
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

    void *donnees;
};

struct AccesseuseAttribut {
    void (*bool_pour_index)(void *, long, bool *);
    void (*entier_pour_index)(void *, long, int *);
    void (*reel_pour_index)(void *, long, float *);
    void (*vec2_pour_index)(void *, long, float *);
    void (*vec3_pour_index)(void *, long, float *);
    void (*vec4_pour_index)(void *, long, float *);
    void (*couleur_pour_index)(void *, long, float *);

    long (*nombre_attributs)(void *);
    void *(*poignee_attribut)(void *, long);
    enum TypeVolume (*type_volume_pour_attribut)(void *);
    void (*nom_attribut)(void *, char **, long *);

    void *donnees_utilisateur;
};

enum MethodeCalculTailleVoxel {
    UNITE_GLOBALE,
    SOUS_DIVISION_AXE_X,
    SOUS_DIVISION_AXE_Y,
    SOUS_DIVISION_AXE_Z,
    SOUS_DIVISION_GRAND_AXE,
};

struct ParametresVDBDepuisMaillage {
    struct AdaptriceMaillage *adaptrice;

    bool utilise_grille_reference;
    struct GrilleVDB *grille_reference;

    bool transfere_attributs;
    bool genere_champs_de_distance;
    bool champs_distance_absolue;
    bool genere_volume_dense;

    struct AccesseuseAttribut *attributs_points;
    struct AccesseuseAttribut *attributs_sommets;
    struct AccesseuseAttribut *attributs_polygones;

    const char *nom_champs_distance;
    const char *nom_volume_dense;

    enum MethodeCalculTailleVoxel methode_calcule_taille_voxel;
    float taille_voxel;
    int compte_voxel;

    bool remplis_interieur;
    bool bande_en_unite_globale;
    float bande_exterieure;
    float bande_interieure;
    int compte_voxel_bande_exterieure;
    int compte_voxel_bande_interieure;
};

struct ContexteEvaluationVDB {
    void (*rapporte_erreur)(void *, const char *, long);
    void (*rapporte_avertissement)(void *, const char *, long);

    void *donnees_utilisateur;
};

void VDB_depuis_polygones(struct ContexteKuri *ctx,
                          struct ContexteEvaluationVDB *ctx_eval,
                          struct ParametresVDBDepuisMaillage *params,
                          struct ExportriceGrilles *exportrice,
                          struct Interruptrice *interruptrice);

#ifdef __cplusplus
}
#endif
