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

#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/ParticlesToLevelSet.h"

using namespace openvdb::OPENVDB_VERSION_NAME;

struct VDBGrid {
};

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
};

#if 0
enum DrapeauxVDBDepuisMaillage {
    CHAMPS_DISTANCE_ABSOLUE = 1,
    DESACTIVE_SUPPRESSION_VOXELS_INTERSECTANT = 2,
    DESACTIVE_RENORMALISATION = 4,
    DESACTIVE_ELAGAGE_LIMITES_BANDE = 8,
};
#endif

struct ParametresVDBDepuisMaillage {
    float taille_voxel;
    float bande_exterieure;
    float bande_interieure;
    bool transfere_attributs;
    bool construit_champs_distance;
    bool champs_distance_absolue;
    bool remplis_interieur_volume;
};

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

VDBGrid *VDB_depuis_polygones(AdaptriceMaillage *adaptrice,
                              ParametresVDBDepuisMaillage *params,
                              Interruptrice *interruptrice)
{
    InterruptriceVDB boss{interruptrice};
    AdaptriceMaillageVDB adaptrice_vdb = AdaptriceMaillageVDB{*adaptrice};
    math::Transform::Ptr transform = math::Transform::createLinearTransform(params->taille_voxel);

    std::vector<Vec3s> pointList = extrait_points(adaptrice_vdb, *transform);
    std::vector<Vec4I> primList = extrait_quads_et_triangles(adaptrice_vdb);

    tools::QuadAndTriangleDataAdapter<Vec3s, Vec4I> mesh(pointList, primList);

    /* Grille d'index pour les attributs. */
    Int32Grid::Ptr grille_index;
    if (params->transfere_attributs) {
        grille_index.reset(new Int32Grid(0));
    }

    int drapeaux = params->champs_distance_absolue ? tools::UNSIGNED_DISTANCE_FIELD : 0;

    FloatGrid::Ptr grille = tools::meshToVolume(boss,
                                                mesh,
                                                *transform,
                                                params->bande_exterieure,
                                                params->bande_interieure,
                                                drapeaux,
                                                grille_index.get());

    if (!boss.wasInterrupted() && params->construit_champs_distance) {
        // À FAIREhvdb::createVdbPrimitive(*gdp, grid, evalStdString("distancename",
        // time).c_str());
    }

    if (!boss.wasInterrupted() && params->remplis_interieur_volume &&
        !params->champs_distance_absolue) {
        // If no level set grid is exported the original level set
        // grid is modified in place.
        FloatGrid::Ptr grille_fog;

        if (params->construit_champs_distance) {
            grille_fog = grille->deepCopy();
        }
        else {
            grille_fog = grille;
        }

        tools::sdfToFogVolume(*grille_fog);

        // À FAIRE hvdb::createVdbPrimitive(*gdp, outputGrid, evalStdString("fogname",
        // time).c_str());
    }

    if (!boss.wasInterrupted() && params->transfere_attributs) {
        // À FAIRE
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

struct AccesseuseAttribut {
    void (*bool_pour_index)(void *, void *, long, bool *);
    void (*entier_pour_index)(void *, void *, long, int *);
    void (*reel_pour_index)(void *, void *, long, float *);
    void (*vec2_pour_index)(void *, void *, long, float *);
    void (*vec3_pour_index)(void *, void *, long, float *);
    void (*vec4_pour_index)(void *, void *, long, float *);
    void (*couleur_pour_index)(void *, void *, long, float *);

    long (*nombre_attributs)(void *);
    void (*nom_attribut)(void *, long, char **, long *);
    void *(*poignee_attribut)(void *, long);
};

struct ParamsVDBDepuisParticules {
    float poids_rayon;
    float poids_velocity;
    bool cree_des_trainees;
    bool transfere_attributs;

    AccesseuseAttribut *acces_attributs_points;
};

VDBGrid *VDB_depuis_particules(AdaptricePoints *adaptrice,
                               ParamsVDBDepuisParticules *params,
                               Interruptrice *interruptrice)
{
    ListeParticules liste_particules = ListeParticules{*adaptrice};
    InterruptriceVDB boss{interruptrice};

    tools::ParticlesToLevelSet p;
}
