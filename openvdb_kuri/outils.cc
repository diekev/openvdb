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

#include "outils.hh"

#include <sstream>

#include "Géométrie3D/ipa.h"

#include "grille_vdb.hh"
#include "ipa_openvdb.h"

namespace outils {

std::string enchaine(std::vector<std::string> const &chaines, std::string const &separateur)
{
    if (chaines.empty()) {
        return "";
    }

    std::stringstream os;

    os << chaines[0];
    for (size_t i = 1; i < chaines.size(); i++) {
        os << separateur;
        os << chaines[i];
    }

    return os.str();
}

std::string chaine_depuis_accesseuse(AccesseuseChaine *accesseuse)
{
    char *nom = nullptr;
    long taille = 0;
    if (accesseuse->accede_chaine) {
        accesseuse->accede_chaine(accesseuse->donnees, &nom, &taille);
    }
    return std::string(nom, static_cast<size_t>(taille));
}

std::vector<GrilleVDB *> grilles_depuis_iteratrice(IteratriceGrillesVDB &iteratrice)
{
    std::vector<GrilleVDB *> resultat;

    auto grille = iteratrice.suivante(iteratrice.donnees_utilisateur);
    while (grille) {
        if (grille->grid) {
            resultat.push_back(grille);
        }
        grille = iteratrice.suivante(iteratrice.donnees_utilisateur);
    }

    return resultat;
}

}  // namespace outils

void InterruptriceVDB::start(const char *name)
{
    if (interruptrice && interruptrice->commence) {
        interruptrice->commence(interruptrice->donnees, name);
    }
}

void InterruptriceVDB::end()
{
    if (interruptrice && interruptrice->termine) {
        interruptrice->termine(interruptrice->donnees);
    }
}

bool InterruptriceVDB::wasInterrupted(int percent)
{
    if (interruptrice && interruptrice->doit_interrompre) {
        return interruptrice->doit_interrompre(interruptrice->donnees, percent);
    }

    return false;
}

size_t AdaptriceMaillageVDB::polygonCount() const
{
    return static_cast<size_t>(this->nombre_de_polygones(this->donnees));
}

size_t AdaptriceMaillageVDB::pointCount() const
{
    return static_cast<size_t>(this->nombre_de_points(this->donnees));
}

size_t AdaptriceMaillageVDB::vertexCount(size_t n) const
{
    return static_cast<size_t>(
        this->nombre_de_sommets_polygone(this->donnees, static_cast<long>(n)));
}

void AdaptriceMaillageVDB::getIndexSpacePoint(size_t n, size_t v, openvdb::math::Vec3d &pos) const
{
    float x, y, z;
    this->point_pour_sommet_polygone(
        this->donnees, static_cast<long>(n), static_cast<long>(v), &x, &y, &z);
    pos.x() = x;
    pos.y() = y;
    pos.z() = z;
}

openvdb::math::Vec3d AdaptriceMaillageVDB::getPoint(long n) const
{
    float x, y, z;
    this->point_pour_index(this->donnees, static_cast<long>(n), &x, &y, &z);
    openvdb::math::Vec3d pos;
    pos.x() = x;
    pos.y() = y;
    pos.z() = z;
    return pos;
}

void AdaptriceMaillageVDB::setPoint(long n, const openvdb::math::Vec3d &point)
{
    this->remplace_point_a_l_index(this->donnees,
                                   n,
                                   static_cast<float>(point.x()),
                                   static_cast<float>(point.y()),
                                   static_cast<float>(point.z()));
}

void AdaptriceMaillageVDB::indexSommetsPolygone(long n, int *index) const
{
    this->index_points_sommets_polygone(this->donnees, n, index);
}

void AdaptriceMaillageVDB::rafinePolygone(long i, const RafineusePolygone &rafineuse) const
{
    if (!this->rafine_polygone) {
        return;
    }
    this->rafine_polygone(this->donnees, i, &rafineuse);
}

openvdb::math::BBox<openvdb::math::Vec3d> AdaptriceMaillageVDB::getBoundBox() const
{
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float min_z = std::numeric_limits<float>::max();
    float max_x = -std::numeric_limits<float>::max();
    float max_y = -std::numeric_limits<float>::max();
    float max_z = -std::numeric_limits<float>::max();
    this->calcule_boite_englobante(this->donnees, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z);
    openvdb::math::Vec3d min = openvdb::math::Vec3d(min_x, min_y, min_z);
    openvdb::math::Vec3d max = openvdb::math::Vec3d(max_x, max_y, max_z);
    return {min, max};
}

openvdb::math::Vec3s AdaptriceMaillageVDB::normalPolygone(size_t i) const
{
    float nx, ny, nz;
    this->calcule_normal_polygone(this->donnees, static_cast<long>(i), &nx, &ny, &nz);
    return {nx, ny, nz};
}

void AdaptriceMaillageVDB::ajoutePoints(float *points, long nombre) const
{
    if (this->ajoute_plusieurs_points) {
        this->ajoute_plusieurs_points(this->donnees, points, nombre);
        return;
    }

    reserveNombreDePoints(nombre);

    for (int i = 0; i < nombre; i++) {
        ajouteUnPoint(points[0], points[1], points[2]);
        points += 3;
    }
}

void AdaptriceMaillageVDB::reserveNombreDePoints(long nombre) const
{
    this->reserve_nombre_de_points(this->donnees, nombre);
}

void AdaptriceMaillageVDB::reserveNombreDePolygones(long nombre) const
{
    this->reserve_nombre_de_polygones(this->donnees, nombre);
}

void AdaptriceMaillageVDB::ajouteUnPoint(float x, float y, float z) const
{
    this->ajoute_un_point(this->donnees, x, y, z);
}

void AdaptriceMaillageVDB::ajouteListePolygones(int *sommets,
                                                int *sommets_par_polygones,
                                                long nombre_polygones)
{
    this->ajoute_liste_polygones(this->donnees, sommets, sommets_par_polygones, nombre_polygones);
}

void AdaptriceMaillageVDB::ajouteUnPolygone(int *sommets, int taille) const
{
    this->ajoute_un_polygone(this->donnees, sommets, taille);
}

void *AdaptriceMaillageVDB::creeUnGroupeDePoints(const std::string &nom) const
{
    return this->cree_un_groupe_de_points(
        this->donnees, nom.c_str(), static_cast<long>(nom.size()));
}

void *AdaptriceMaillageVDB::creeUnGroupeDePolygones(const std::string &nom) const
{
    return this->cree_un_groupe_de_polygones(
        this->donnees, nom.c_str(), static_cast<long>(nom.size()));
}

void AdaptriceMaillageVDB::ajouteAuGroupe(void *poignee_groupe, long index) const
{
    this->ajoute_au_groupe(poignee_groupe, index);
}

void AdaptriceMaillageVDB::ajoutePlageAuGroupe(void *poignee_groupe,
                                               long index_debut,
                                               long index_fin) const
{
    this->ajoute_plage_au_groupe(poignee_groupe, index_debut, index_fin);
}

bool AdaptriceMaillageVDB::groupePolygonePossedePoint(const void *poignee_groupe, long index) const
{
    return this->groupe_polygone_possede_point(poignee_groupe, index);
}

tbb::blocked_range<long> AdaptriceMaillageVDB::plagePoint() const
{
    return tbb::blocked_range<long>(0, static_cast<long>(pointCount()));
}
