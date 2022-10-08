/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2022 Kévin Dietrich. */

#pragma once

#include "Géométrie3D/ipa.h"

#include <string>
#include <vector>

#include <openvdb/math/BBox.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointDataGrid.h>

#include <tbb/blocked_range.h>

struct AccesseuseChaine;
struct Interruptrice;
struct ContexteKuri;
struct ExportriceGrilles;
struct GrilleVDB;
struct IteratriceGrillesVDB;

namespace outils {

std::string enchaine(std::vector<std::string> const &chaines, std::string const &separateur);

std::string chaine_depuis_accesseuse(AccesseuseChaine *accesseuse);

std::vector<GrilleVDB *> grilles_depuis_iteratrice(IteratriceGrillesVDB &iteratrice);

void exporte_grille_vdb(ContexteKuri *ctx,
                        ExportriceGrilles *exportrice,
                        openvdb::GridBase::Ptr grille_vdb,
                        std::string nom);

}  // namespace outils

/* *********************************************************************** */

// Grid type lists, for use with GEO_PrimVDB::apply(), GEOvdbApply(),
// or openvdb::GridBase::apply()

using ScalarGridTypes = openvdb::TypeList<openvdb::BoolGrid,
                                          openvdb::FloatGrid,
                                          openvdb::DoubleGrid,
                                          openvdb::Int32Grid,
                                          openvdb::Int64Grid>;

using PointGridTypes = openvdb::TypeList<openvdb::points::PointDataGrid>;

using VolumeGridTypes = ScalarGridTypes::Append<openvdb::Vec3GridTypes>;

using AllGridTypes = VolumeGridTypes::Append<PointGridTypes>;

/* *********************************************************************** */

/* Enveloppe autour de l'Interruptrice afin de la conformer à l'interface attendue par OpenVDB. */
struct InterruptriceVDB {
    Interruptrice *interruptrice = nullptr;

    void start(const char *name = nullptr);

    void end();

    bool wasInterrupted(int percent = -1);
};

/* *********************************************************************** */

using PlageDonnees = tbb::blocked_range<long>;

inline bool est_valide(const PlageDonnees plage)
{
    return plage.begin() <= plage.end();
}

// interface pour OpenVDB
class AdaptriceMaillageVDB : public AdaptriceMaillage {
  public:
    static AdaptriceMaillageVDB enveloppe(AdaptriceMaillage *adaptrice)
    {
        return {*adaptrice};
    }

    size_t polygonCount() const;

    size_t pointCount() const;

    size_t vertexCount(size_t n) const;

    void getIndexSpacePoint(size_t n, size_t v, openvdb::math::Vec3d &pos) const;

    openvdb::math::Vec3d getPoint(long n) const;

    void setPoint(long n, openvdb::math::Vec3d const &point);

    void indexSommetsPolygone(long n, int *index) const;

    void rafinePolygone(long i, const RafineusePolygone &rafineuse) const;

    openvdb::math::BBox<openvdb::math::Vec3d> getBoundBox() const;

    openvdb::math::Vec3s normalPolygone(size_t i) const;

    void ajoutePoints(float *points, long nombre) const;

    void reserveNombreDePoints(long nombre) const;

    void reserveNombreDePolygones(long nombre) const;

    void ajouteUnPoint(float x, float y, float z) const;

    void ajouteListePolygones(int *sommets, int *sommets_par_polygones, long nombre_polygones);

    void ajouteUnPolygone(int *sommets, int taille) const;

    void *creeUnGroupeDePoints(const std::string &nom) const;

    void *creeUnGroupeDePolygones(const std::string &nom) const;

    void ajouteAuGroupe(void *poignee_groupe, long index) const;

    void ajoutePlageAuGroupe(void *poignee_groupe, long index_debut, long index_fin) const;

    bool groupePolygonePossedePoint(const void *poignee_groupe, long index) const;

    PlageDonnees plagePoint(long decalage = 0) const;
};

/* *********************************************************************** */

struct EnveloppeContexteEvaluation : public ContexteEvaluation {
    static EnveloppeContexteEvaluation enveloppe(ContexteEvaluation *ctx)
    {
        return {*ctx};
    }

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
