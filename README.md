# omniCADD
#### a CADD scoring system to asses variant deleteriousness in non-model organisms
----

## Overview
This repo contains the reconstructed workflow of computing CADD (Combined Annotation Dependent Depletion) scores, as first described in [Kircher et al 2014](https://www.nature.com/articles/ng.2892). 

The code it an update, modification and extension based on the CADD models for [mouse](https://doi.org/10.1186/s12859-018-2337-5) (Groß et al 2018) , [pig](https://doi.org/10.1186/s12711-020-0528-9) (Groß et al 2020a) and [chicken](https://doi.org/10.1371/journal.pgen.1009027) (Groß et al 2020b)

### main pipeline - dense data
As of now, the main branch contains the workflow for running it using already deposited data, such as EPO alignments from Ensembl, deposited reference genomes, and Ensembl VEP. It also assumed that the user has a population-wide VCF file with variants to be used in the section *derive variants*. 

### alignment free pipeline - scarce data
The *alignmentfree* branch (still under development), rather, assumes the opposite, namely that

- The species of interest is not part of a publicly available alignment
- the reference species is not deposited in Ensembl
- The species of interest does not have annotation (gff\[3\])

The plan is to merge and make the pipeline fully modifiable in the future. 
#### :exclamation: **NOTE** :exclamation: 
The scripts are still under construction and have not yet fully been wrapped into the workflow. 

## :link: user prerequisities 
To work correctly, the pipeline assumes that critical files, such as the reference genome, epo alignment and population-level-VCF are in subfolders in `resources`. Before running, make sure they are either transfered there, or downloaded directly there before running. 
```bash
cd resources/genome
# download here, example workflow has been run with domesticated pig:
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz

cd ../alignment
# download here, example workflow has been run with the epo 44 way mammal alignment:
wget -r -nd --no-parent -e robots=off -A '44_mammals.epo.*.maf.gz' https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/44_mammals.epo/

cd ../pop-level-VCF
# transfer here, example workflow has been run with local pig vcf files
```

## :ballot_box_with_check: TO-DO

- clean up scripts and workflow
- rempove temporary information
- make cluster profile

## Information
*TBA*

## Usage
*TBA*

## References
Kircher, M., Witten, D., Jain, P. et al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet 46, 310–315 (2014). doi.org/10.1038/ng.2892

Groß, C., de Ridder, D. & Reinders, M. Predicting variant deleteriousness in non-human species: applying the CADD approach in mouse. BMC Bioinformatics 19, 373 (2018). doi.org/10.1186/s12859-018-2337-5

Groß, C., Derks, M., Megens, HJ. et al. pCADD: SNV prioritisation in Sus scrofa. Genet Sel Evol 52, 4 (2020). doi.org/10.1186/s12711-020-0528-9

Groß C, Bortoluzzi C, de Ridder D, Megens HJ, Groenen MAM, et al. (2020) Prioritizing sequence variants in conserved non-coding elements in the chicken genome using chCADD. PLOS Genetics 16(9): e1009027. doi.org/10.1371/journal.pgen.1009027

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor.
Genome Biology Jun 6;17(1):122. (2016)
doi:10.1186/s13059-016-0974-4 