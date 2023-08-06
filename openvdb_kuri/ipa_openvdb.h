/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2022 Kévin Dietrich. */

#pragma once

#ifdef __cplusplus
extern "C" {
#else
typedef char bool;
#endif

struct AdaptriceMaillage;
struct ContexteEvaluation;
struct ContexteKuri;
struct GrilleVDB;
struct Interruptrice;

enum TypeVolume {
    VOLUME_INVALIDE = 0,
    VOLUME_R32 = 1,
    VOLUME_R64 = 2,
    VOLUME_Z32 = 3,
    VOLUME_Z64 = 4,
    VOLUME_BOOL = 5,
    VOLUME_VEC3_R32 = 6,
    VOLUME_VEC3_R64 = 7,
    VOLUME_VEC3_Z32 = 8,
    VOLUME_INDEX_POINT = 9,
    VOLUME_DONNEES_POINT = 10,
};

struct GrilleVDB *VDB_copie_grille(struct ContexteKuri *ctx, struct GrilleVDB *grille);
void VDB_detruit_grille(struct ContexteKuri *ctx, struct GrilleVDB *grille);
void VDB_donne_nom_grille(struct GrilleVDB *grille, const char **nom, long *taille);
void VDB_definis_nom_grille(struct GrilleVDB *grille, const char *nom, long taille);
enum TypeVolume VDB_type_volume_pour_grille(struct GrilleVDB *grille);

struct ExportriceGrilles {
    void (*ajoute_grille)(void *, struct GrilleVDB *);
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

struct AccesseuseChaine {
    void (*accede_chaine)(void *, char **, long *);
    void *donnees;
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

    struct AccesseuseChaine *nom_champs_distance;
    struct AccesseuseChaine *nom_volume_dense;

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

void VDB_depuis_polygones(struct ContexteKuri *ctx,
                          struct ContexteEvaluation *ctx_eval,
                          struct ParametresVDBDepuisMaillage *params,
                          struct ExportriceGrilles *exportrice,
                          struct Interruptrice *interruptrice);

struct IteratriceGrillesVDB {
    struct GrilleVDB *(*suivante)(void *);

    void *donnees_utilisateur;
};

struct ParametresVDBVersMaillage {
    struct AdaptriceMaillage *maillage_reference;
    struct IteratriceGrillesVDB *groupe_grilles;

    struct GrilleVDB *grille_masque_surface;
    struct GrilleVDB *grille_champs_adaptivite;

    float isovalue;
    float adaptivite;
    float decalage_masque;
    bool inverse_masque;

    /* Options pour la surface de référence. */

    /**
     * Quand une surface de référence est fournie, ceci définie le seuil d'adaptivité pour les
     * régions qui se trouvent à l'interieur de la surface (par exemple, les faces internes
     * générées par une fracture de la surface d'origine).
     * @min 0.0
     * @max 1.0
     */
    float adaptivite_interne;

    /** Affine les bords et coins du maillage généré.
     * @défaut vrai
     */
    bool affiner_les_traits;

    /** Controle le masque d'adaptivité des bords.
     * @min 0.0
     * @max 1.0
     * @défaut 0.5
     */
    float tolerance_de_bord;
};

void VDB_vers_polygones(struct ContexteKuri *ctx,
                        struct ContexteEvaluation *ctx_eval,
                        struct ParametresVDBVersMaillage *params,
                        struct AdaptriceMaillage *maillage,
                        struct Interruptrice *interruptrice);

enum ComportementFichierManquant {
    RAPPORTE_ERREUR,
    CONTINUE_SANS_RIEN_CREER,
};

struct ParametresLectureVDB {
    /** Chemin vers le fichier `.vdb`. */
    struct AccesseuseChaine *chemin_fichier;

    /** Détermine que faire si le fichier n'est pas trouvé : émettre une erreur et arrêter
     * (RAPPORTE_ERREUR), ou continue sans rien créer (CONTINUE_SANS_RIEN_CREER). */
    enum ComportementFichierManquant comportement_fichier_manquant;

    /** Crée des grilles ne contenant que les métadonnées des grilles dans le fichier. */
    bool metadonnees_seules;

    /** Rogne les grilles selon la boite englobante du maillage de référence. */
    bool rogne;

    /** Maillage de référence pour déterminer la boite englobante. */
    struct AdaptriceMaillage *maillage_reference;

    /** Crée un groupe contenant les grilles lues. */
    bool cree_groupe;

    /** Nom du groupe à créer, si #cree_groupe est actif. */
    struct AccesseuseChaine *nom_groupe;

    /** Ne charge les données des grilles que lorsqu'elles seront vraiment utilisées. */
    bool chargement_tardif;

    /** Si #chargement_tardif est actif, les fichiers plus petit que cette limite (en giga-octets)
     * sont copier dans un dossier temporaire afin de garantir que le fichier ne sera pas modifié
     * par un programme externe. */
    float limite_pour_copier;

    // À FAIRE : liste de grilles à charger
};

void VDB_depuis_fichier(struct ContexteKuri *ctx,
                        struct ContexteEvaluation *ctx_eval,
                        struct ParametresLectureVDB *params,
                        struct ExportriceGrilles *flux_sortie_grille,
                        struct Interruptrice *interruptrice);

/* ------------------------------------------------------------------------- */
/** \name Écriture de fichier .vdb
 * \{ */

enum MethodeCompressionVDB {
    METHODE_COMPRESSION_VDB_AUCUNE,
    METHODE_COMPRESSION_VDB_ZIP,
    METHODE_COMPRESSION_VDB_BLOSC,
};

struct ParametresEcritureVDB {
    /** Chemin vers le fichier de sortie. */
    struct AccesseuseChaine *nom_fichier_sortie;

    /** Chaine pour définir un attribut dans le fichier ayant le nom de l'application ayant créée
     * ledit fichier. */
    struct AccesseuseChaine *nom_application;

    /** Ensemble des grilles à écrire. */
    struct IteratriceGrillesVDB *grilles;

    /** Définis la méthode de compression. ZIP est lent mais comprime très bien. Blosc est rapide
     * et comprime bien. Pour la plupart des cas, Blosc est recommandé. */
    enum MethodeCompressionVDB methode_compression;

    /** Pour les grilles appartenant à ce groupe, les valeurs scalaires ou vectorielles en point
     * flottant seront écrites en 16-bit. Si aucune grille n'est présente dans ce groupe, les
     * grilles seront écrites selon leurs propres paramètres de précision. Les grilles ici doivent
     * aussi se trouver dans les grilles à écrire. */
    struct IteratriceGrillesVDB *grilles_precision_16_bit;

    /** Pour les grilles appartenant à ce groupe, les valeurs scalaires ou vectorielles en point
     * flottant seront écrites en 32-bit. Si aucune grille n'est présente dans ce groupe, les
     * grilles seront écrites selon leurs propres paramètres de précision. Les grilles ici doivent
     * aussi se trouver dans les grilles à écrire. */
    struct IteratriceGrillesVDB *grilles_precision_32_bit;
};

void VDB_ecris_fichier(struct ContexteEvaluation *ctx_eval,
                       struct ParametresEcritureVDB *params,
                       struct Interruptrice *interruptrice);

/** \} */

/* *************************************************************************** */
/* Distribution de points.
 */

enum PolitiqueNommageGrille {
    GARDE_NOM_ORIGINAL,
    AJOUTE_SUFFIXE,
    REMPLACE_NOM,
};

enum DomaineChampsDeDistance {
    INTERIEUR,
    SURFACE,
    BANDE_FINE,
};

enum CompressionPoints {
    NE_COMPRIME_PAS,
    COMPRIME_SUR_16_BIT,
    COMPRIME_SUR_8_BIT,
};

enum ModeGenereationPoints {
    PAR_COMPTE_TOTAL,
    PAR_DENSITE_LOCALE,
    PAR_VOXEL,
};

struct ParametresPointsDepuisVDB {
    struct AdaptriceMaillage *sortie_points;
    struct IteratriceGrillesVDB *groupe_grilles;
    struct ExportriceGrilles *sortie_grilles;

    /**
     * Génère une grille de points OpenVDB au lieu de points natifs.
     */
    bool genere_grille_points;

    bool rogne_selon_surface;

    bool verbeux;

    int graine;

    float diffusion;

    int nombre_de_points;

    float points_par_voxel;

    float densite;

    bool multiplie_par_densite;
    enum PolitiqueNommageGrille nommage_grilles_sortie;

    enum DomaineChampsDeDistance domaine_champs_de_distance;

    struct AccesseuseChaine *nom_grille_sortie;

    enum CompressionPoints compression_desiree;

    enum ModeGenereationPoints mode_generation_points;

    float valeur_iso;

    bool cree_groupe;

    struct AccesseuseChaine *nom_groupe_sortie;
};

void VDB_distribue_points(struct ContexteKuri *ctx,
                          struct ContexteEvaluation *ctx_eval,
                          struct ParametresPointsDepuisVDB *params,
                          struct Interruptrice *interruptrice);

struct ParamsVDBDepuisPoints {
    float poids_rayon;
    float poids_velocity;
    bool cree_des_trainees;
    bool transfere_attributs;

    struct AccesseuseAttribut *acces_attributs_points;
};

void VDB_depuis_points(struct ContexteKuri *ctx,
                       struct ContexteEvaluation *ctx_eval,
                       struct AdaptricePoints *adaptrice,
                       struct ParamsVDBDepuisPoints *params,
                       struct Interruptrice *interruptrice);

/* ------------------------------------------------------------------------- */
/** \name Dimension mondiale d'une grille.
 * \{ */

struct DimensionMondeVDB {
    float min_x;
    float min_y;
    float min_z;

    float max_x;
    float max_y;
    float max_z;
};

void VDB_donne_dimension_monde(struct GrilleVDB *grille, struct DimensionMondeVDB *r_dimension);

/** \} */

/* ------------------------------------------------------------------------- */
/** \name Conversion vers tampon dense.
 * \{ */

struct DonneesConversionVersDense {
    float *donnees;
    int taille_donnees;
    int dim_x;
    int dim_y;
    int dim_z;
};

void VDB_donne_tampon_dense(struct ContexteKuri *ctx,
                            struct GrilleVDB *grille,
                            struct DonneesConversionVersDense *r_donnees);

void VDB_detruit_tampon_dense(struct ContexteKuri *ctx,
                              struct DonneesConversionVersDense *donnees);

/** \} */

#ifdef __cplusplus
}
#endif
