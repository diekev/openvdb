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
};

VDBGrid *VDB_depuis_polygones(AdaptriceMaillage *adaptrice)
{
    AdaptriceMaillageVDB adaptrice_vdb = AdaptriceMaillageVDB{*adaptrice};
    math::Transform transform;
    float exterior_bandwidth = 3.0f;
    float interiorBandWidth = 3.0f;
    int flags = 0;

    tools::meshToVolume(adaptrice_vdb, transform, exterior_bandwidth, interiorBandWidth, flags);
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

VDBGrid *VDB_depuis_particules(AdaptricePoints *adaptrice,
                               ParamsVDBDepuisParticules *params,
                               Interruptrice *interruptrice)
{
    ListeParticules liste_particules = ListeParticules{*adaptrice};
    InterruptriceVDB boss{interruptrice};

    tools::ParticlesToLevelSet p;
}
