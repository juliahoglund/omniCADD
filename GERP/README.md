# How to cactus
#### https://github.com/harvardinformatics/GenomeAnnotation-WholeGenomeAlignment

0. load modules needed to run it in a container:
```bash
ml bioinfo-tools mafTools
ml PDC apptainer hal
module unload python
ml python3
# ml cactus
```

1. pull the (CPU) image, unless GPUs are available
```bash
aptainer pull --disable-cache docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.0
```

2. chekc tree and what genomic fastas are needed,
then download (here: [43 amniotes](https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/43_mammals.epo/README.43_mammals.epo))
(i am sure there is a better way, but for now, single commands.)
```bash
wget https://ftp.ensembl.org/pub/current_fasta/loxodonta_africana/dna/Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/marmota_marmota_marmota/dna/Marmota_marmota_marmota.marMar2.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/sciurus_vulgaris/dna/Sciurus_vulgaris.mSciVul1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cricetulus_griseus_chok1gshd/dna/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/microtus_ochrogaster/dna/Microtus_ochrogaster.MicOch1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/peromyscus_maniculatus_bairdii/dna/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_spretus/dna/Mus_spretus.SPRET_EiJ_v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_caroli/dna/Mus_caroli.CAROLI_EIJ_v1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_pahari/dna/Mus_pahari.PAHARI_EIJ_v1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cavia_porcellus/dna/Cavia_porcellus.Cavpor3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/microcebus_murinus/dna/Microcebus_murinus.Mmur_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/pan_paniscus/dna/Pan_paniscus.panpan1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/gorilla_gorilla/dna/Gorilla_gorilla.gorGor4.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/nomascus_leucogenys/dna/Nomascus_leucogenys.Nleu_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/canis_lupus_dingo/dna/Canis_lupus_dingo.ASM325472v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/panthera_leo/dna/Panthera_leo.PanLeo1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/panthera_pardus/dna/Panthera_pardus.PanPar1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/felis_catus/dna/Felis_catus.Felis_catus_9.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_indicus_hybrid/dna/Bos_indicus_hybrid.UOA_Brahman_1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_grunniens/dna/Bos_grunniens.LU_Bosgru_v3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/capra_hircus/dna/Capra_hircus.ARS1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/ovis_aries_rambouillet/dna/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cervus_hanglu_yarkandensis/dna/Cervus_hanglu_yarkandensis.CEY_v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/physeter_catodon/dna/Physeter_catodon.ASM283717v2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/phocoena_sinus/dna/Phocoena_sinus.mPhoSin1.pri.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/delphinapterus_leucas/dna/Delphinapterus_leucas.ASM228892v3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/monodon_monoceros/dna/Monodon_monoceros.NGI_Narwhal_1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/balaenoptera_musculus/dna/Balaenoptera_musculus.mBalMus1.v2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/catagonus_wagneri/dna/Catagonus_wagneri.CatWag_v2_BIUU_UCD.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/camelus_dromedarius/dna/Camelus_dromedarius.CamDro2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/rhinolophus_ferrumequinum/dna/Rhinolophus_ferrumequinum.mRhiFer1_v1.p.dna_sm.toplevel.fa.gz


```

3. decide how new genome is added, and create the corresponding info files for cactus.
here, asian elephan (elephas maximus) is being added to the same node as the afrixan elephant 
(loxofdonta africana) and will probably be the outermost specied, i.e. there will be only one subtree
corresponding to the old alignment, or maybe african eephant as the only leaf in the top alignment.
-----

add-to-node.txt
```
(((((((((cricetulus_griseus_chok1gshd:0.06086,microtus_ochrogaster:0.06521)Cricetidae:0.00471,peromyscus_maniculatus_bairdii:0.05821)Cricetidae:0.01299,((((mus_spretus:0.01485,mus_musculus:0.01891)Mus:0.00973,mus_caroli:0.02376)Mus:0.01384,mus_pahari:0.03623)Mus:0.02324,rattus_norvegicus:0.05417)Murinae:0.03551)Muroidea:0.08994,cavia_porcellus:0.14738)Rodentia:0.00405,(sciurus_vulgaris:0.04795,marmota_marmota_marmota:0.05572)Sciuridae:0.04079)Rodentia:0.01251,oryctolagus_cuniculus:0.13568)Glires:0.00663,((((((pan_troglodytes:0.00712,pan_paniscus:0.00937)Pan:0.00585,homo_sapiens:0.01021)Homininae:0.00331,gorilla_gorilla:0.01986)Homininae:0.01348,nomascus_leucogenys:0.03718)Hominoidea:0.00852,((macaca_mulatta:0.00903,macaca_fascicularis:0.01613)Macaca:0.01584,chlorocebus_sabaeus:0.01721)Cercopithecinae:0.02009)Catarrhini:0.05414,microcebus_murinus:0.08543)Primates:0.00973)Euarchontoglires:0.012,((((((panthera_pardus:0.00778,panthera_leo:0.00922)Panthera:0.0099,felis_catus:0.01636)Felidae:0.04768,(canis_lupus_dingo:0.00521,canis_lupus_familiaris:0.00596)Canis_lupus:0.06758)Carnivora:0.02451,equus_caballus:0.06742)Laurasiatheria:0.00386,(((((((monodon_monoceros:0.00782,delphinapterus_leucas:0.00873)Monodontidae:0.00695,phocoena_sinus:0.01691)Odontoceti:0.01694,physeter_catodon:0.03414)Odontoceti:0.00425,balaenoptera_musculus:0.0269)Cetacea:0.02935,(((ovis_aries_rambouillet:0.01416,capra_hircus:0.01621)Caprinae:0.01847,((bos_taurus:0.00588,bos_indicus_hybrid:0.00675)Bos:0.0044,bos_grunniens:0.02771)Bos:0.0201)Bovidae:0.00995,cervus_hanglu_yarkandensis:0.03491)Pecora:0.04656)Artiodactyla:0.01025,(sus_scrofa:0.04408,catagonus_wagneri:0.04764)Suina:0.03591)Artiodactyla:0.00614,camelus_dromedarius:0.07127)Artiodactyla:0.01949)Laurasiatheria:0.00321,rhinolophus_ferrumequinum:0.09309)Laurasiatheria:0.01407)Boreoeutheria:0.00982,(elephas_maximus,loxodonta_africana:0.11174)Proboscidae)Eutheria;

cricetulus_griseus_chok1gshd	Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa.gz
microtus_ochrogaster	Microtus_ochrogaster.MicOch1.0.dna_sm.toplevel.fa.gz
peromyscus_maniculatus_bairdii	Peromyscus_maniculatus_bairdii.HU_Pman_2.1.dna_sm.toplevel.fa.gz
mus_spretus	Mus_spretus.SPRET_EiJ_v1.dna_sm.toplevel.fa.gz
mus_musculus	Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz
mus_caroli	Mus_caroli.CAROLI_EIJ_v1.1.dna_sm.toplevel.fa.gz
mus_pahari	Mus_pahari.PAHARI_EIJ_v1.1.dna_sm.toplevel.fa.gz
rattus_norvegicus	Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.gz
cavia_porcellus	Cavia_porcellus.Cavpor3.0.dna_sm.toplevel.fa.gz
sciurus_vulgaris	Sciurus_vulgaris.mSciVul1.1.dna_sm.toplevel.fa.gz
marmota_marmota_marmota	Marmota_marmota_marmota.marMar2.1.dna_sm.toplevel.fa.gz
oryctolagus_cuniculus	Oryctolagus_cuniculus.OryCun2.0.dna_sm.toplevel.fa.gz
pan_troglodytes	Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa.gz
pan_paniscus	Pan_paniscus.panpan1.1.dna_sm.toplevel.fa.gz
homo_sapiens	Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
gorilla_gorilla	Gorilla_gorilla.gorGor4.dna_sm.toplevel.fa.gz
nomascus_leucogenys	Nomascus_leucogenys.Nleu_3.0.dna_sm.toplevel.fa.gz
macaca_mulatta	Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.gz
macaca_fascicularis	Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz
chlorocebus_sabaeus	Chlorocebus_sabaeus.ChlSab1.1.dna_sm.toplevel.fa.gz
microcebus_murinus	Microcebus_murinus.Mmur_3.0.dna_sm.toplevel.fa.gz
panthera_pardus	Panthera_pardus.PanPar1.0.dna_sm.toplevel.fa.gz
panthera_leo	Panthera_leo.PanLeo1.0.dna_sm.toplevel.fa.gz
felis_catus	Felis_catus.Felis_catus_9.0.dna_sm.toplevel.fa.gz
canis_lupus_dingo	Canis_lupus_dingo.ASM325472v1.dna_sm.toplevel.fa.gz
canis_lupus_familiaris	Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz
equus_caballus	Equus_caballus.EquCab3.0.dna_sm.toplevel.fa.gz
monodon_monoceros	Monodon_monoceros.NGI_Narwhal_1.dna_sm.toplevel.fa.gz
delphinapterus_leucas	Delphinapterus_leucas.ASM228892v3.dna_sm.toplevel.fa.gz
phocoena_sinus	Phocoena_sinus.mPhoSin1.pri.dna_sm.toplevel.fa.gz
physeter_catodon	Physeter_catodon.ASM283717v2.dna_sm.toplevel.fa.gz
balaenoptera_musculus	Balaenoptera_musculus.mBalMus1.v2.dna_sm.toplevel.fa.gz
ovis_aries_rambouillet	Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_sm.toplevel.fa.gz
capra_hircus	Capra_hircus.ARS1.dna_sm.toplevel.fa.gz
bos_taurus	Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz
bos_grunniens	Bos_grunniens.LU_Bosgru_v3.0.dna_sm.toplevel.fa.gz
cervus_hanglu_yarkandensis	Cervus_hanglu_yarkandensis.CEY_v1.dna_sm.toplevel.fa.gz
sus_scrofa	Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz
catagonus_wagneri	Catagonus_wagneri.CatWag_v2_BIUU_UCD.dna_sm.toplevel.fa.gz
camelus_dromedarius	Camelus_dromedarius.CamDro2.dna_sm.toplevel.fa.gz
rhinolophus_ferrumequinum	Rhinolophus_ferrumequinum.mRhiFer1_v1.p.dna_sm.toplevel.fa.gz
loxodonta_africana	Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa.gz
elephas_maximus	GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna.gz

```

4. do the actual preparation to add a node
```bash
cactus-update-prepare add node ./evolverMammals.hal ./input.txt --genome mr --outDir ./steps --jobStore ./jobstore
cactus-update-prepare add node [hal alignment] input-seq.txt --genome Eutheria --outDir ./steps --jobStore ./jobstore

# Specifically, for the Singularity image, it can be run as:

singularity exec --cleanenv cactus_v2.2.0-gpu.sif cactus-prepare <INPUT FILE> --outDir <OUTPUT DIRECTORY> --jobStore <TEMP DIRECTORY> --gpu

## to run
singularity exec --bind results/alignment/hal/43_mammals.epo:/mnt cactus_v2.9.0.sif ls /mnt
singularity exec --cleanenv cactus_v2.9.0.sif cactus-update-prepare add node $alignment input-seq.txt --genome Eutheria --outDir ./steps --jobStore ./jobstore

# and then add to make it work
singularity exec --cleanenv --bind $PWD:/mnt cactus_v2.9.0.sif cactus-update-prepare add node /mnt/results/alignment/hal/43_mammals.epo/43_mammals.epo.1_1.hal /mnt/input-seq.txt --genome Eutheria --outDir /mnt/steps/ --jobStore /mnt/jobstore/

halStats results/alignment/hal/43_mammals.epo/43_mammals.epo.1_1.hal

# OMG HOLY FUCK	!!!
https://ftp.ensembl.org/pub/misc/compara/multi/hal_files/Mammals-100-way_20230606.hal
# 253 gb /\
https://ftp.ensembl.org/pub/current_compara/species_trees/100_mammals_Cactus_default.nh

# singularity exec --cleanenv --bind $PWD:/mnt cactus_v2.9.0.sif cactus-update-prepare add branch --parentGenome loxodonta_africana --childGenome cricetulus_griseus_chok1gshd /mnt/results/alignment/hal/43_mammals.epo/43_mammals.epo.1_1.hal /mnt/input-seq.txt --cactus-prepare-options '--alignCores 4' --outDir /mnt/steps/ --jobStore /mnt/jobstore/ --ancestorName AncElephant 
```

