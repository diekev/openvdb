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

#include "ipa_openvdb.h"

#include "openvdb/tools/LevelSetUtil.h"
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/ParticlesToLevelSet.h"

#include "../InterfaceCKuri/contexte_kuri.hh"

using namespace openvdb::OPENVDB_VERSION_NAME;

/* *********************************************************************** */

struct GrilleVDB {
    GridBase::Ptr grid;
};

void VDB_detruit_grille(struct ContexteKuri *ctx, struct GrilleVDB *grille)
{
    kuri_deloge(ctx, grille);
}

void VDB_accede_nom_grille(struct GrilleVDB *grille, const char **nom, long *taille)
{
    if (!grille || !grille->grid) {
        *nom = nullptr;
        *taille = 0;
        return;
    }

    const auto &nom_grille = grille->grid->getName();
    *nom = nom_grille.c_str();
    *taille = static_cast<long>(nom_grille.size());
}

void VDB_mute_nom_grille(struct GrilleVDB *grille, const char *nom, long taille)
{
    if (!grille || !grille->grid) {
        return;
    }

    const auto nouveau_nom = std::string(nom, static_cast<size_t>(taille));
    grille->grid->setName(nouveau_nom);
}

TypeVolume VDB_type_volume_pour_grille(GrilleVDB *grille)
{
    if (!grille) {
        return TypeVolume::INVALIDE;
    }

    GridBase::ConstPtr grille_vdb = grille->grid;
    if (!grille_vdb) {
        return TypeVolume::INVALIDE;
    }

#define VERIFIE_TYPE_ET_RETOURNE(GridType, type_volume)                                           \
    if (grille_vdb->isType<GridType>()) {                                                         \
        return TypeVolume::type_volume;                                                           \
    }

    VERIFIE_TYPE_ET_RETOURNE(FloatGrid, R32)
    VERIFIE_TYPE_ET_RETOURNE(DoubleGrid, R64)
    VERIFIE_TYPE_ET_RETOURNE(Int32Grid, Z32)
    VERIFIE_TYPE_ET_RETOURNE(Int64Grid, Z64)
    VERIFIE_TYPE_ET_RETOURNE(BoolGrid, BOOL)
    VERIFIE_TYPE_ET_RETOURNE(Vec3SGrid, VEC3_R32)
    VERIFIE_TYPE_ET_RETOURNE(Vec3DGrid, VEC3_R64)
    VERIFIE_TYPE_ET_RETOURNE(Vec3IGrid, VEC3_Z32)
    VERIFIE_TYPE_ET_RETOURNE(tools::PointIndexGrid, INDEX_POINT)
    VERIFIE_TYPE_ET_RETOURNE(points::PointDataGrid, DONNEES_POINT)

#undef VERIFIE_TYPE_ET_RETOURNE

    return TypeVolume::INVALIDE;
}

/* *********************************************************************** */

static void exporte_grille_vdb(ContexteKuri *ctx,
                               ExportriceGrilles *exportrice,
                               GridBase::Ptr grille_vdb,
                               std::string nom)
{
    grille_vdb->setName(nom);
    GrilleVDB *grille = kuri_loge<GrilleVDB>(ctx);
    grille->grid = grille_vdb;
    exportrice->ajoute_grille(exportrice->donnees, grille);
}

/* *********************************************************************** */

struct InterruptriceVDB {
    Interruptrice *interruptrice = nullptr;

    void start(const char *name = nullptr)
    {
        if (interruptrice && interruptrice->commence) {
            interruptrice->commence(interruptrice->donnees, name);
        }
    }

    void end()
    {
        if (interruptrice && interruptrice->termine) {
            interruptrice->termine(interruptrice->donnees);
        }
    }

    bool wasInterrupted(int percent = -1)
    {
        if (interruptrice && interruptrice->doit_interrompre) {
            return interruptrice->doit_interrompre(interruptrice->donnees, percent);
        }

        return false;
    }
};

// interface pour OpenVDB
class AdaptriceMaillageVDB : public AdaptriceMaillage {
  public:
    size_t polygonCount() const
    {
        return static_cast<size_t>(this->nombre_de_polygones(this->donnees));
    }

    size_t pointCount() const
    {
        return static_cast<size_t>(this->nombre_de_points(this->donnees));
    }

    size_t vertexCount(size_t n) const
    {
        return static_cast<size_t>(
            this->nombre_de_sommets_polygone(this->donnees, static_cast<long>(n)));
    }

    void getIndexSpacePoint(size_t n, size_t v, math::Vec3d &pos) const
    {
        float x, y, z;
        this->point_pour_sommet_polygone(
            this->donnees, static_cast<long>(n), static_cast<long>(v), &x, &y, &z);
        pos.x() = x;
        pos.y() = y;
        pos.z() = z;
    }

    math::Vec3d getPoint(long n) const
    {
        float x, y, z;
        this->point_pour_index(this->donnees, static_cast<long>(n), &x, &y, &z);
        math::Vec3d pos;
        pos.x() = x;
        pos.y() = y;
        pos.z() = z;
        return pos;
    }

    void indexSommetsPolygone(long n, int *index) const
    {
        this->index_points_sommets_polygone(this->donnees, n, index);
    }

    void rafinePolygone(long i, const RafineusePolygone &rafineuse) const
    {
        this->rafine_polygone(this->donnees, i, &rafineuse);
    }

    math::BBox<Vec3d> getBoundBox() const
    {
        float min_x = std::numeric_limits<float>::max();
        float min_y = std::numeric_limits<float>::max();
        float min_z = std::numeric_limits<float>::max();
        float max_x = -std::numeric_limits<float>::max();
        float max_y = -std::numeric_limits<float>::max();
        float max_z = -std::numeric_limits<float>::max();
        this->calcule_boite_englobante(
            this->donnees, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z);
        Vec3d min = Vec3d(min_x, min_y, min_z);
        Vec3d max = Vec3d(max_x, max_y, max_z);
        return {min, max};
    }
};

static AdaptriceMaillageVDB enveloppe(AdaptriceMaillage *adaptrice)
{
    return {*adaptrice};
}

#if 0
enum DrapeauxVDBDepuisMaillage {
    CHAMPS_DISTANCE_ABSOLUE = 1,
    DESACTIVE_SUPPRESSION_VOXELS_INTERSECTANT = 2,
    DESACTIVE_RENORMALISATION = 4,
    DESACTIVE_ELAGAGE_LIMITES_BANDE = 8,
};
#endif

static float calcule_taille_voxel_effective(const ParametresVDBDepuisMaillage &params)
{
    auto taille_voxel_depuis_compte = [&](float taille_axe) {
        const float dim = static_cast<float>(params.compte_voxel);
        if (params.bande_en_unite_globale) {
            const float w = params.bande_exterieure;
            return (taille_axe + 2.0f * w) / dim;
        }

        const float w = static_cast<float>(params.compte_voxel_bande_interieure);
        return taille_axe / std::max(1.0f, dim - 2.0f * w);
    };

    switch (params.methode_calcule_taille_voxel) {
        case MethodeCalculTailleVoxel::UNITE_GLOBALE:
        {
            return params.taille_voxel;
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_X:
        {
            const auto bbox = enveloppe(params.adaptrice).getBoundBox();
            const float taille_x = bbox.extents().x();
            return taille_voxel_depuis_compte(taille_x);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_Y:
        {
            const auto bbox = enveloppe(params.adaptrice).getBoundBox();
            const float taille_y = bbox.extents().y();
            return taille_voxel_depuis_compte(taille_y);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_AXE_Z:
        {
            const auto bbox = enveloppe(params.adaptrice).getBoundBox();
            const float taille_z = bbox.extents().z();
            return taille_voxel_depuis_compte(taille_z);
        }
        case MethodeCalculTailleVoxel::SOUS_DIVISION_GRAND_AXE:
        {
            const auto bbox = enveloppe(params.adaptrice).getBoundBox();
            const float taille_max = std::max(bbox.extents().x(),
                                              std::max(bbox.extents().y(), bbox.extents().z()));
            return taille_voxel_depuis_compte(taille_max);
        }
    }

    return 0.0f;
}

static float calcule_taille_bande_effective(const ParametresVDBDepuisMaillage &params,
                                            float bande_globale,
                                            float bande_voxel,
                                            float taille_voxel)
{
    if (params.bande_en_unite_globale) {
        return bande_globale / taille_voxel;
    }
    return bande_voxel;
}

struct DonneesCreationGrilleAttribut {
    std::string nom_attribut;
    void *poignee_attribut;
    TypeVolume type_volume;
};

static std::vector<DonneesCreationGrilleAttribut> parse_attributs_a_transferer(
    AccesseuseAttribut *accesseuse)
{
    std::vector<DonneesCreationGrilleAttribut> resultat;

    const long nombre_attributs = accesseuse->nombre_attributs(accesseuse->donnees_utilisateur);
    if (nombre_attributs == 0) {
        return resultat;
    }

    for (int i = 0; i < nombre_attributs; i++) {
        void *poignee = accesseuse->poignee_attribut(accesseuse->donnees_utilisateur, i);
        if (!poignee) {
            continue;
        }

        TypeVolume type_volume = accesseuse->type_volume_pour_attribut(poignee);
        if (type_volume == TypeVolume::INVALIDE) {
            continue;
        }

        char *nom_attribut = nullptr;
        long taille_nom_attribut = 0;
        accesseuse->nom_attribut(poignee, &nom_attribut, &taille_nom_attribut);
    }

    return resultat;
}

static std::vector<Vec3s> extrait_points(const AdaptriceMaillageVDB &adaptrice_vdb,
                                         const math::Transform &transform)
{
    std::vector<Vec3s> resultat;
    const long nombre_de_points = adaptrice_vdb.pointCount();

    if (nombre_de_points == 0) {
        return resultat;
    }

    resultat.reserve(nombre_de_points);

    for (int i = 0; i < nombre_de_points; i++) {
        Vec3d point = adaptrice_vdb.getPoint(i);
        point = transform.worldToIndex(point);
        resultat.push_back(point);
    }

    return resultat;
}

static std::vector<Vec4I> extrait_quads_et_triangles(const AdaptriceMaillageVDB &adaptrice_vdb)
{
    std::vector<Vec4I> resultat;
    const long nombre_de_primitives = adaptrice_vdb.polygonCount();

    if (nombre_de_primitives == 0) {
        return resultat;
    }

    resultat.reserve(nombre_de_primitives);

    RafineusePolygone rafineuse;
    rafineuse.donnees = &resultat;
    rafineuse.ajoute_triangle = [](RafineusePolygone *raf, long v0, long v1, long v2) {
        Vec4I triangle;
        triangle.x() = static_cast<int>(v0);
        triangle.y() = static_cast<int>(v1);
        triangle.z() = static_cast<int>(v2);
        triangle.w() = util::INVALID_IDX;
        static_cast<std::vector<Vec4I> *>(raf->donnees)->push_back(triangle);
    };
    rafineuse.ajoute_quadrilatere =
        [](RafineusePolygone *raf, long v0, long v1, long v2, long v3) {
            Vec4I triangle;
            triangle.x() = static_cast<int>(v0);
            triangle.y() = static_cast<int>(v1);
            triangle.z() = static_cast<int>(v2);
            triangle.w() = static_cast<int>(v3);
            static_cast<std::vector<Vec4I> *>(raf->donnees)->push_back(triangle);
        };

    for (int i = 0; i < nombre_de_primitives; i++) {
        const long nombre_de_sommets = adaptrice_vdb.vertexCount(i);

        if (nombre_de_sommets == 3) {
            int index[3];
            adaptrice_vdb.indexSommetsPolygone(i, index);
            Vec4I triangle;
            triangle.x() = index[0];
            triangle.y() = index[1];
            triangle.z() = index[2];
            triangle.w() = util::INVALID_IDX;
            resultat.push_back(triangle);
        }
        else if (nombre_de_sommets == 4) {
            int index[4];
            adaptrice_vdb.indexSommetsPolygone(i, index);
            Vec4I quad;
            quad.x() = index[0];
            quad.y() = index[1];
            quad.z() = index[2];
            quad.w() = index[3];
            resultat.push_back(quad);
        }
        else if (nombre_de_sommets > 4) {
            if (!adaptrice_vdb.rafine_polygone) {
                continue;
            }

            adaptrice_vdb.rafinePolygone(i, rafineuse);
        }
    }

    return resultat;
}

struct EnveloppeContexteEvaluationVDB : public ContexteEvaluationVDB {
    template <typename... Args>
    void rapporteErreur(Args... args)
    {
        std::ostringstream ostr;
        ((ostr << args), ...);
        this->rapporteErreur(ostr.str());
    }

    void rapporteErreur(const std::string &erreur)
    {
        if (!this->rapporte_erreur) {
            return;
        }

        this->rapporte_erreur(
            this->donnees_utilisateur, erreur.c_str(), static_cast<long>(erreur.size()));
    }
    template <typename... Args>
    void rapporteAvertissement(Args... args)
    {
        std::ostringstream ostr;
        ((ostr << args), ...);
        this->rapporteAvertissement(ostr.str());
    }

    void rapporteAvertissement(const std::string &avertissement)
    {
        if (!this->rapporte_avertissement) {
            return;
        }
        this->rapporte_avertissement(this->donnees_utilisateur,
                                     avertissement.c_str(),
                                     static_cast<long>(avertissement.size()));
    }
};

static EnveloppeContexteEvaluationVDB enveloppe(ContexteEvaluationVDB *ctx)
{
    return {*ctx};
}

static void VDB_depuis_polygones_impl(ContexteKuri *ctx,
                                      EnveloppeContexteEvaluationVDB &ctx_eval_,
                                      ParametresVDBDepuisMaillage *params,
                                      ExportriceGrilles *exportrice,
                                      InterruptriceVDB &boss)
{
    boss.start("Conversion de polygones en volumes VDB");

    if (!params->genere_champs_de_distance && !params->genere_volume_dense &&
        !params->transfere_attributs) {
        ctx_eval_.rapporteErreur("Rien à générer");
        return;
    }

    math::Transform::Ptr transform;
    float taille_voxel;

    if (params->utilise_grille_reference) {
        if (!params->grille_reference || !params->grille_reference->grid) {
            ctx_eval_.rapporteErreur("Aucune grille de reference trouvée");
            return;
        }

        transform = params->grille_reference->grid->transform().copy();
        taille_voxel = static_cast<float>(transform->voxelSize()[0]);
    }
    else {
        taille_voxel = calcule_taille_voxel_effective(*params);
        transform = math::Transform::createLinearTransform(taille_voxel);
    }

    if (taille_voxel < 1e-5) {
        ctx_eval_.rapporteErreur("La taille de voxel (", taille_voxel, ") est trop petite");
        return;
    }

    AdaptriceMaillageVDB adaptrice_vdb = enveloppe(params->adaptrice);

    std::vector<Vec3s> points = extrait_points(adaptrice_vdb, *transform);
    std::vector<Vec4I> polygones = extrait_quads_et_triangles(adaptrice_vdb);

    /* Même si nous n'avons aucun point ou polygone, nous créons les grilles, car des opérations
     * suivantes pourraient tout de même gérer des grilles vides. */
    if (points.empty()) {
        ctx_eval_.rapporteAvertissement("Aucun point dans le maillage d'entrée");
    }
    if (polygones.empty()) {
        ctx_eval_.rapporteAvertissement("Aucun polygone à convertir");
    }

    tools::QuadAndTriangleDataAdapter<Vec3s, Vec4I> mesh(points, polygones);

    /* Grille d'index pour les attributs. */
    Int32Grid::Ptr grille_index;
    if (params->transfere_attributs) {
        grille_index.reset(new Int32Grid(0));
    }

    int drapeaux = params->champs_distance_absolue ? tools::UNSIGNED_DISTANCE_FIELD : 0;

    const float bande_exterieure = calcule_taille_bande_effective(
        *params, params->bande_exterieure, params->compte_voxel_bande_exterieure, taille_voxel);
    const float bande_interieure = params->remplis_interieur ?
                                       std::numeric_limits<float>::max() :
                                       calcule_taille_bande_effective(
                                           *params,
                                           params->bande_interieure,
                                           params->compte_voxel_bande_interieure,
                                           taille_voxel);

    FloatGrid::Ptr grille = tools::meshToVolume<FloatGrid>(
        boss, mesh, *transform, bande_exterieure, bande_interieure, drapeaux, grille_index.get());

    if (!boss.wasInterrupted() && params->genere_champs_de_distance) {
        exporte_grille_vdb(ctx, exportrice, grille, params->nom_champs_distance);
    }

    if (params->genere_volume_dense && params->champs_distance_absolue) {
        ctx_eval_.rapporteAvertissement("Impossible de générer un volume dense car le champs de "
                                        "distance à générer possède des valeurs absolues");
    }

    if (!boss.wasInterrupted() && params->genere_volume_dense &&
        !params->champs_distance_absolue) {
        /* Si nous n'exportons pas le champs de distance, nous pouvons le modifier directemet. */
        FloatGrid::Ptr grille_fog;
        if (params->genere_champs_de_distance) {
            grille_fog = grille->deepCopy();
        }
        else {
            grille_fog = grille;
        }

        tools::sdfToFogVolume(*grille_fog);

        exporte_grille_vdb(ctx, exportrice, grille_fog, params->nom_volume_dense);
    }

    if (!boss.wasInterrupted() && params->transfere_attributs) {
        // À FAIRE
    }

    if (boss.wasInterrupted()) {
        ctx_eval_.rapporteAvertissement("Le processus fut interrompu");
    }

    boss.end();
}

void VDB_depuis_polygones(ContexteKuri *ctx,
                          ContexteEvaluationVDB *ctx_eval,
                          ParametresVDBDepuisMaillage *params,
                          ExportriceGrilles *exportrice,
                          Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    EnveloppeContexteEvaluationVDB ctx_eval_ = enveloppe(ctx_eval);

    try {
        VDB_depuis_polygones_impl(ctx, ctx_eval_, params, exportrice, boss);
    }
    catch (std::exception &e) {
        ctx_eval_.rapporteErreur(e.what());
    }
}

struct ListeParticules : public AdaptricePoints {
    // Interface de convenience.
    bool possede_rayon() const
    {
        return this->possede_attribut_pour_rayon &&
               this->possede_attribut_pour_rayon(this->donnees);
    }

    bool possede_velocite() const
    {
        return this->possede_attribut_pour_velocite &&
               this->possede_attribut_pour_velocite(this->donnees);
    }

    // Interface requise par OpenVDB.
    size_t size() const
    {
        return static_cast<size_t>(this->nombre_de_points(this->donnees));
    }

    void getPos(size_t n, Vec3R &xyz) const
    {
        float x, y, z;
        this->position_pour_index(this->donnees, static_cast<size_t>(n), &x, &y, &z);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
    }

    void getPosRad(size_t n, Vec3R &xyz, Real &radius) const
    {
        float x, y, z, rayon;
        this->position_et_rayon_pour_index(
            this->donnees, static_cast<size_t>(n), &x, &y, &z, &rayon);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
        radius = rayon;
    }

    void getPosRadVel(size_t n, Vec3R &xyz, Real &radius, Vec3R &velocity) const
    {
        float x, y, z, rayon, vx, vy, vz;
        this->position_rayon_et_velocite_pour_index(
            this->donnees, static_cast<size_t>(n), &x, &y, &z, &rayon, &vx, &vy, &vz);
        xyz.x() = x;
        xyz.y() = y;
        xyz.z() = z;
        radius = rayon;
        velocity.x() = vx;
        velocity.y() = vy;
        velocity.z() = vz;
    }

    /* Pour le transfert d'attributs, nous utilisons une grille d'index qui sera ensuite utilisée
     * pour transférer les attributs. Cette grille stocke l'index du point le plus proche de chaque
     * voxel. */
    void getAtt(size_t n, Int32 &att) const
    {
        att = static_cast<Int32>(n);
    }
};

struct ParamsVDBDepuisParticules {
    float poids_rayon;
    float poids_velocity;
    bool cree_des_trainees;
    bool transfere_attributs;

    AccesseuseAttribut *acces_attributs_points;
};

void VDB_depuis_particules(AdaptricePoints *adaptrice,
                           ParamsVDBDepuisParticules *params,
                           Interruptrice *interruptrice)
{
    //    ListeParticules liste_particules = ListeParticules{*adaptrice};
    //    InterruptriceVDB boss{interruptrice};

    //    tools::ParticlesToLevelSet p;
}
