/* SPDX-License-Identifier: GPL-2.0-or-later
 * The Original Code is Copyright (C) 2023 Kévin Dietrich. */

#include "lecture_vdb.hh"

#include "grille_vdb.hh"
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

static std::set<GrilleVDB const *> ensemble_grilles_depuis_iteratrice(IteratriceGrillesVDB *iter)
{
    if (!iter) {
        return {};
    }

    std::set<GrilleVDB const *> résultat;

    auto grilles = outils::grilles_depuis_iteratrice(*iter);
    for (auto grille : grilles) {
        résultat.insert(grille);
    }

    return résultat;
}

static bool ensemble_possède(std::set<GrilleVDB const *> const &ensemble, GrilleVDB const *grille)
{
    return ensemble.find(grille) != ensemble.end();
}

static uint32_t détermine_compression_fichier(io::File const &fichier,
                                              MethodeCompressionVDB méthode_compression)
{

    uint32_t compression = fichier.compression();
    switch (méthode_compression) {
        case METHODE_COMPRESSION_VDB_AUCUNE:
        {
            compression &= ~(io::COMPRESS_ZIP | io::COMPRESS_BLOSC);
            break;
        }
        case METHODE_COMPRESSION_VDB_ZIP:
        {
#ifdef OPENVDB_USE_ZLIB
            compression |= io::COMPRESS_ZIP;
            compression &= ~io::COMPRESS_BLOSC;
#else
            compression = io::COMPRESS_ACTIVE_MASK;
#endif
            break;
        }
        case METHODE_COMPRESSION_VDB_BLOSC:
        {
#ifdef OPENVDB_USE_BLOSC
            compression &= ~io::COMPRESS_ZIP;
            compression |= io::COMPRESS_BLOSC;
#else
            compression = io::COMPRESS_ACTIVE_MASK;
#endif
            break;
        }
    }

    return compression;
}

static void rapporte_conflits_de_précision(EnveloppeContexteEvaluation &ctx_eval,
                                           std::set<std::string> const &conflits_précision)
{
    if (conflits_précision.empty()) {
        return;
    }

    std::stringstream ss;
    if (conflits_précision.size() == 1) {
        ss << "La grille '" << *conflits_précision.begin() << "' se trouve ";
    }
    else {
        ss << "Les grilles ";

        std::string virgule = "'";

        for (auto const &nom : conflits_précision) {
            ss << virgule << nom << "'";
            virgule = ", '";
        }

        ss << " se trouvent ";
    }

    ss << "à la fois dans le groupe de précision 16-bit et dans celui de 32-bit. Choisissez-en un "
          "mais pas les deux.";

    ctx_eval.rapporteAvertissement(ss.str());
}

void écriture_vdb(EnveloppeContexteEvaluation &ctx_eval,
                  ParametresEcritureVDB *params,
                  InterruptriceVDB &boss)
{
    auto chemin_sortie = outils::chaine_depuis_accesseuse(params->nom_fichier_sortie);
    if (chemin_sortie.empty()) {
        ctx_eval.rapporteErreur("Aucun fichier de sortie précisé");
        return;
    }

    if (!params->grilles) {
        ctx_eval.rapporteAvertissement("Aucune grille trouvée");
        return;
    }

    auto grilles_à_écrire = outils::grilles_depuis_iteratrice(*params->grilles);
    if (grilles_à_écrire.empty()) {
        ctx_eval.rapporteAvertissement("Aucune grille trouvée");
        return;
    }

    auto grilles_precision_16_bit = ensemble_grilles_depuis_iteratrice(
        params->grilles_precision_16_bit);
    auto grilles_precision_32_bit = ensemble_grilles_depuis_iteratrice(
        params->grilles_precision_32_bit);

    std::set<std::string> conflits_précision;

    GridPtrSet grilles_sorties;
    for (auto grille : grilles_à_écrire) {
        if (boss.wasInterrupted()) {
            return;
        }

        if (!grille->grid) {
            continue;
        }

        /* Copie pour pouvoir changer les paramètres d'écriture entre autre. */
        auto grid = grille->grid->copyGrid();

        // À FAIRE : nommage des grilles.
        auto const nom_grille = grid->getName();

        auto const écris_en_16_bit = ensemble_possède(grilles_precision_16_bit, grille);
        auto const écris_en_32_bit = ensemble_possède(grilles_precision_32_bit, grille);

        if (écris_en_16_bit && écris_en_32_bit) {
            conflits_précision.insert(nom_grille);
        }
        else if (écris_en_16_bit) {
            grid->setSaveFloatAsHalf(true);
        }
        else if (écris_en_32_bit) {
            grid->setSaveFloatAsHalf(false);
        }

        grilles_sorties.insert(grid);
    }

    rapporte_conflits_de_précision(ctx_eval, conflits_précision);

    /* Ajout des métadonnées du fichier. */
    MetaMap métadonnées_fichier;
    auto nom_application = outils::chaine_depuis_accesseuse(params->nom_application);
    if (nom_application != "") {
        métadonnées_fichier.insertMeta("creator", StringMetadata(nom_application));
    }

    /* Création du fichier. */
    io::File fichier(chemin_sortie);

    auto const compression = détermine_compression_fichier(fichier, params->methode_compression);
    fichier.setCompression(compression);

    fichier.write(grilles_sorties, métadonnées_fichier);
    fichier.close();
}
}  // namespace kvdb
