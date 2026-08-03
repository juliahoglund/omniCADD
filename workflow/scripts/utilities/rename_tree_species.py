#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename the leaf labels of the 44-mammal Ensembl species tree from Ensembl's
naming convention (Genus_species_AssemblyOrStrainName) to the lowercase
genus_species-style names used for the alignment's species labels
throughout this pipeline (e.g. clean_ambiguous, mark_ancestor, sort_by_chr).

Without this, gerpcol/phyloFit/phastCons see zero overlap between the
alignment's species labels and the tree's leaf labels, since the raw
Ensembl download uses a completely different naming convention.

The mapping below is not a generic pattern - several species need
different, inconsistent amounts of their name kept (e.g.
Marmota_marmota_marmota_marMar2.1 keeps all three repeated words while
Gorilla_gorilla_gorilla_gorGor4 drops to two), so it is hardcoded and
verified against the actual alignment species list rather than derived
from a regex.

Usage: python3 rename_tree_species.py INPUT.nwk OUTPUT.nwk
"""

import sys

# tree leaf label -> alignment species label
RENAME = {
    "Cricetulus_griseus_CHOK1GS_HDv1": "cricetulus_griseus_chok1gshd",
    "Microtus_ochrogaster_MicOch1.0": "microtus_ochrogaster",
    "Peromyscus_maniculatus_bairdii_HU_Pman_2.1": "peromyscus_maniculatus_bairdii",
    "Mus_musculus_reference_CL57BL6_strain": "mus_musculus",
    "Mus_caroli_strain_CAROLI_EIJ": "mus_caroli",
    "Mus_pahari_strain_PAHARI_EIJ": "mus_pahari",
    "Rattus_norvegicus_reference_strain": "rattus_norvegicus",
    "Cavia_porcellus_Cavpor3.0": "cavia_porcellus",
    "Marmota_marmota_marmota_marMar2.1": "marmota_marmota_marmota",
    "Sciurus_vulgaris_mSciVul1.1": "sciurus_vulgaris",
    "Oryctolagus_cuniculus_OryCun2.0": "oryctolagus_cuniculus",
    "Pan_troglodytes_Pan_tro_3.0": "pan_troglodytes",
    "Pan_paniscus_panpan1.1": "pan_paniscus",
    "Homo_sapiens_GRCh38": "homo_sapiens",
    "Gorilla_gorilla_gorilla_gorGor4": "gorilla_gorilla",
    "Pongo_abelii_Susie_PABv2": "pongo_abelii",
    "Nomascus_leucogenys_Nleu_3.0": "nomascus_leucogenys",
    "Macaca_mulatta_Mmul_10": "macaca_mulatta",
    "Macaca_fascicularis_Macaca_fascicularis_6.0": "macaca_fascicularis",
    "Papio_anubis_Panubis1.0": "papio_anubis",
    "Chlorocebus_sabaeus_ChlSab1.1": "chlorocebus_sabaeus",
    "Callithrix_jacchus_mCalJac1.pat.X": "callithrix_jacchus",
    "Microcebus_murinus_Mmur_3.0": "microcebus_murinus",
    "Monodon_monoceros_NGI_Narwhal_1": "monodon_monoceros",
    "Delphinapterus_leucas_ASM228892v3": "delphinapterus_leucas",
    "Phocoena_sinus_mPhoSin1.pri": "phocoena_sinus",
    # historical synonym: Physeter catodon / Physeter macrocephalus (sperm whale)
    "Physeter_macrocephalus_ASM283717v2": "physeter_catodon",
    "Balaenoptera_musculus_mBalMus1.v2": "balaenoptera_musculus",
    "Ovis_aries_reference_breed": "ovis_aries",
    "Capra_hircus_reference_breed": "capra_hircus",
    "Bos_taurus_reference_breed": "bos_taurus",
    "Bos_indicus_x_Bos_taurus_UOA_Brahman_1": "bos_indicus_hybrid",
    "Bos_grunniens_LU_Bosgru_v3.0": "bos_grunniens",
    "Cervus_hanglu_yarkandensis_CEY_v1": "cervus_hanglu_yarkandensis",
    "Catagonus_wagneri_CatWag_v2_BIUU_UCD": "catagonus_wagneri",
    "Sus_scrofa_reference_breed": "sus_scrofa",
    "Camelus_dromedarius_CamDro2": "camelus_dromedarius",
    "Rhinolophus_ferrumequinum_mRhiFer1_v1.p": "rhinolophus_ferrumequinum",
    "Panthera_pardus_PanPar1.0": "panthera_pardus",
    "Panthera_leo_PanLeo1.0": "panthera_leo",
    "Canis_lupus_familiaris_reference_breed": "canis_lupus_familiaris",
    "Canis_lupus_dingo_ASM325472v1": "canis_lupus_dingo",
    "Equus_caballus_breed_thoroughbred": "equus_caballus",
    "Loxodonta_africana_loxAfr3": "loxodonta_africana",
}


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: rename_tree_species.py INPUT.nwk OUTPUT.nwk")

    input_path, output_path = sys.argv[1], sys.argv[2]

    with open(input_path, "r") as f:
        nwk = f.read()

    missing = [leaf for leaf in RENAME if leaf not in nwk]
    if missing:
        sys.exit(f"ERROR: expected leaf label(s) not found in {input_path}: {missing}")

    for leaf, alignment_name in RENAME.items():
        nwk = nwk.replace(leaf, alignment_name)

    with open(output_path, "w") as f:
        f.write(nwk)

    print(f"Renamed {len(RENAME)} leaf labels, wrote {output_path}")


if __name__ == "__main__":
    main()
