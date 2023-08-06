/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "vdb_depuis_points.hh"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"
#include "outils.hh"

using namespace openvdb;

namespace kvdb {

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

void vdb_depuis_points(ContexteKuri *ctx,
                       EnveloppeContexteEvaluation &ctx_eval,
                       AdaptricePoints *adaptrice,
                       ParamsVDBDepuisPoints *params,
                       InterruptriceVDB &interruptrice)
{
    //    ListeParticules liste_particules = ListeParticules{*adaptrice};

    //    tools::ParticlesToLevelSet p;
}

}  // namespace kvdb
