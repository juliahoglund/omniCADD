#!/bin/bash
#
# Download mammalian genomes for multiway alignment generation
# Uses Ensembl FTP to download 43 species from EPO alignment
#
# Author: Julia Höglund
# Date: 2026-02-27
#

set -euo pipefail

# Configuration
GENOME_DIR="resources/alignment_genomes"
BASE_URL="https://ftp.ensembl.org/pub/current_fasta"
LOG_FILE="${GENOME_DIR}/download.log"

# Create output directory
mkdir -p "$GENOME_DIR"

# Start logging
echo "Starting genome downloads at $(date)" | tee -a "$LOG_FILE"
echo "Output directory: $GENOME_DIR" | tee -a "$LOG_FILE"

# Function to download genome
download_genome() {
    local species=$1
    local assembly=$2
    local url="${BASE_URL}/${species}/dna/${species}.${assembly}.dna_sm.toplevel.fa.gz"
    local output="${GENOME_DIR}/${species}.fa.gz"
    
    if [ -f "$output" ]; then
        echo "SKIP: $species (already downloaded)" | tee -a "$LOG_FILE"
        return 0
    fi
    
    echo "Downloading: $species ($assembly)" | tee -a "$LOG_FILE"
    if wget -q --show-progress -O "$output" "$url" 2>> "$LOG_FILE"; then
        echo "SUCCESS: $species" | tee -a "$LOG_FILE"
    else
        echo "FAILED: $species" | tee -a "$LOG_FILE"
        rm -f "$output"
        return 1
    fi
}

# Download all species (43 mammals from EPO alignment)
echo "Downloading 43 mammalian genomes..." | tee -a "$LOG_FILE"

# Afrotheria
download_genome "loxodonta_africana" "loxAfr3"

# Lagomorpha
download_genome "oryctolagus_cuniculus" "OryCun2.0"

# Rodentia (Sciuromorpha)
download_genome "marmota_marmota_marmota" "marMar2.1"
download_genome "sciurus_vulgaris" "mSciVul1.1"

# Rodentia (Myomorpha)
download_genome "cricetulus_griseus_chok1gshd" "CHOK1GS_HDv1"
download_genome "microtus_ochrogaster" "MicOch1.0"
download_genome "peromyscus_maniculatus_bairdii" "HU_Pman_2.1"
download_genome "mus_musculus" "GRCm39"
download_genome "mus_spretus" "SPRET_EiJ_v1"
download_genome "mus_caroli" "CAROLI_EIJ_v1.1"
download_genome "mus_pahari" "PAHARI_EIJ_v1.1"
download_genome "rattus_norvegicus" "mRatBN7.2"

# Rodentia (Hystricomorpha)
download_genome "cavia_porcellus" "Cavpor3.0"

# Primates (Strepsirrhini)
download_genome "microcebus_murinus" "Mmur_3.0"

# Primates (Catarrhini)
download_genome "chlorocebus_sabaeus" "ChlSab1.1"
download_genome "macaca_mulatta" "Mmul_10"
download_genome "macaca_fascicularis" "Macaca_fascicularis_6.0"
download_genome "homo_sapiens" "GRCh38"
download_genome "pan_paniscus" "panpan1.1"
download_genome "pan_troglodytes" "Pan_tro_3.0"
download_genome "gorilla_gorilla" "gorGor4"
download_genome "nomascus_leucogenys" "Nleu_3.0"

# Carnivora (Caniformia)
download_genome "canis_lupus_dingo" "ASM325472v1"
download_genome "canis_lupus_familiaris" "ROS_Cfam_1.0"

# Carnivora (Feliformia)
download_genome "panthera_leo" "PanLeo1.0"
download_genome "panthera_pardus" "PanPar1.0"
download_genome "felis_catus" "Felis_catus_9.0"

# Perissodactyla
download_genome "equus_caballus" "EquCab3.0"

# Artiodactyla (Ruminantia)
download_genome "bos_indicus_hybrid" "UOA_Brahman_1"
download_genome "bos_taurus" "ARS-UCD1.3"
download_genome "bos_grunniens" "LU_Bosgru_v3.0"
download_genome "capra_hircus" "ARS1"
download_genome "ovis_aries_rambouillet" "ARS-UI_Ramb_v2.0"
download_genome "cervus_hanglu_yarkandensis" "CEY_v1"

# Cetacea
download_genome "physeter_catodon" "ASM283717v2"
download_genome "phocoena_sinus" "mPhoSin1.pri"
download_genome "delphinapterus_leucas" "ASM228892v3"
download_genome "monodon_monoceros" "NGI_Narwhal_1"
download_genome "balaenoptera_musculus" "mBalMus1.v2"

# Artiodactyla (Suina)
download_genome "sus_scrofa" "Sscrofa11.1"
download_genome "catagonus_wagneri" "CatWag_v2_BIUU_UCD"

# Tylopoda
download_genome "camelus_dromedarius" "CamDro2"

# Chiroptera
download_genome "rhinolophus_ferrumequinum" "mRhiFer1_v1.p"

# Summary
echo "" | tee -a "$LOG_FILE"
echo "Download complete at $(date)" | tee -a "$LOG_FILE"
echo "Downloaded genomes: $(ls -1 ${GENOME_DIR}/*.fa.gz 2>/dev/null | wc -l)" | tee -a "$LOG_FILE"
echo "Total size: $(du -sh ${GENOME_DIR} | cut -f1)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Next steps:" | tee -a "$LOG_FILE"
echo "  1. Verify downloads: ls -lh ${GENOME_DIR}" | tee -a "$LOG_FILE"
echo "  2. Run alignment workflow:" | tee -a "$LOG_FILE"
echo "     snakemake -s workflow/run_alignment_generation.smk \\" | tee -a "$LOG_FILE"
echo "         --configfile config/alignment_generation.yaml \\" | tee -a "$LOG_FILE"
echo "         --cores 32 --use-conda" | tee -a "$LOG_FILE"
