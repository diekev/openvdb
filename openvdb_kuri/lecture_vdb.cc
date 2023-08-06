/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "lecture_vdb.hh"

#include "ipa_openvdb.h"
#include "outils.hh"

using namespace openvdb;

namespace kvdb {

void lecture_vdb(ContexteKuri *ctx,
                 EnveloppeContexteEvaluation &ctx_eval,
                 ParametresLectureVDB *params,
                 ExportriceGrilles *flux_sortie_grille,
                 InterruptriceVDB &boss)
{
    boss.start("Lecture d'un fichier OpenVDB");

    const std::string chemin_fichier = outils::chaine_depuis_accesseuse(params->chemin_fichier);
    const Index64 limite_copie = static_cast<Index64>(1.0e9 * params->limite_pour_copier);

    io::File fichier(chemin_fichier);
    fichier.setCopyMaxBytes(limite_copie);
    fichier.open(params->chargement_tardif);

    MetaMap::Ptr metadonnees_fichier = fichier.getMetadata();
    if (!metadonnees_fichier) {
        metadonnees_fichier.reset(new MetaMap);
    }

    bool rogne = params->rogne;
    BBoxd boite_rognage;
    if (rogne && params->maillage_reference) {
        AdaptriceMaillageVDB maillage = AdaptriceMaillageVDB::enveloppe(
            params->maillage_reference);
        boite_rognage = maillage.getBoundBox();
        rogne = boite_rognage.isSorted();
    }

    for (io::File::NameIterator iter_nom = fichier.beginName(); iter_nom != fichier.endName();
         ++iter_nom) {
        if (boss.wasInterrupted()) {
            throw std::runtime_error("Le processus fut interrompu");
        }

        const std::string &nom_grille = iter_nom.gridName();

        GridBase::Ptr grille;
        if (params->metadonnees_seules) {
            grille = fichier.readGridMetadata(nom_grille);
        }
        else if (rogne) {
            grille = fichier.readGrid(nom_grille, boite_rognage);
        }
        else {
            grille = fichier.readGrid(nom_grille);
        }

        if (!grille) {
            continue;
        }

        /* Transfère les métadonnées du fichier à la grille. */
        for (MetaMap::ConstMetaIterator iter_metadonnee = metadonnees_fichier->beginMeta(),
                                        end = metadonnees_fichier->endMeta();
             iter_metadonnee != end;
             ++iter_metadonnee) {
            /* Ne transfère une métadonnée que si la grille n'en a pas une avec un nom similaire.
             */
            if (Metadata::Ptr meta = iter_metadonnee->second) {
                const std::string name = iter_metadonnee->first;
                if (!(*grille)[name]) {
                    grille->insertMeta(name, *meta);
                }
            }
        }

        outils::exporte_grille_vdb(ctx, flux_sortie_grille, grille, nom_grille);
    }
    fichier.close();

    boss.end();
}
}  // namespace kvdb
