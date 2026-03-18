omniCADD Documentation
======================

**omniCADD** is a CADD scoring system to assess variant deleteriousness in non-model organisms.

.. image:: https://readthedocs.org/projects/omnicadd/badge/?version=latest
    :target: https://omnicadd.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Overview
--------

omniCADD is a **Snakemake-based bioinformatics pipeline** that computes CADD (Combined Annotation Dependent Depletion) 
scores—a measure of variant deleteriousness—for **non-model organisms**.

The pipeline combines:

* **Evolutionary information** from multi-species alignments
* **Functional annotations** (VEP/SNPEff)
* **Conservation metrics** (PHAST, GERP)
* **Machine learning** (logistic regression)

To generate **PHRED-scaled deleteriousness scores** for every possible variant in a genome.

This workflow is a reconstruction of the CADD scoring system, as first described in 
`Kircher et al 2014 <https://www.nature.com/articles/ng.2892>`_.

The code is an update, modification and extension based on the CADD models for:

* `Mouse <https://doi.org/10.1186/s12859-018-2337-5>`_ (Groß et al 2018)
* `Pig <https://doi.org/10.1186/s12711-020-0528-9>`_ (Groß et al 2020a)
* `Chicken <https://doi.org/10.1371/journal.pgen.1009027>`_ (Groß et al 2020b)

Quick Links
-----------

* **GitHub Repository**: https://github.com/juliahoglund/omniCADD
* **Issue Tracker**: https://github.com/juliahoglund/omniCADD/issues
* **Documentation**: You're reading it!

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   
   installation
   
.. toctree::
   :maxdepth: 2
   :caption: Setup Guides
   
   SETUP
   DARDEL_SETUP
   DOCKER_SLURM_GUIDE

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   workflow
   configuration
   input_data
   output_data
   
.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics
   
   DOCKER_SLURM_GUIDE
   
.. toctree::
   :maxdepth: 1
   :caption: Reference
   
   workflow

Workflow Overview
-----------------

omniCADD processes genomic data through 8 main steps:

1. **Extract Ancestral Sequences** - Reconstruct ancestral sequences from multi-species alignments
2. **Derive Variants** - Identify real variants from population VCF data
3. **Simulate Variants** - Generate neutral variants for model training
4. **Summary Report** - Create interactive HTML visualization of genome and variants
5. **Annotate Variants** - Apply functional (VEP) and conservation (PHAST, GERP) annotations
6. **Combine Annotations** - Merge all annotations into feature matrices
7. **Train & Test Model** - Build logistic regression classifier with cross-validation
8. **Score Variants** - Generate genome-wide PHRED-scaled CADD scores

Pipeline Modes
--------------

Main Pipeline - Dense Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main branch processes species with rich genomic resources:

* EPO alignments from Ensembl
* Deposited reference genomes
* Ensembl VEP annotation
* Population-wide VCF files with allele frequencies

This mode provides the highest quality CADD scores.

Alignment-Free Pipeline - Scarce Data  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *alignmentfree* branch (under development) supports species with limited resources:

* No public multi-species alignment
* Reference not in Ensembl
* Limited or no gene annotation

**Status**: Under development, not yet fully integrated

Requirements
------------

Software Dependencies
~~~~~~~~~~~~~~~~~~~~~

* **Snakemake** ≥ 7.0
* **Conda/Mamba** for environment management
* **Docker** or **Singularity/Apptainer** (optional, for containerized execution)

System Requirements
~~~~~~~~~~~~~~~~~~~

* **RAM**: 16+ GB recommended (64+ GB for whole genomes)
* **Storage**: 500+ GB for alignments, annotations, and intermediate files
* **CPU**: Multi-core processor recommended for parallelization

Input Data
----------

Required Inputs
~~~~~~~~~~~~~~~

The pipeline requires the following data in the ``resources/`` directory:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1
   
   * - Data Type
     - Location
     - Description
   * - Reference Genome
     - ``resources/genome/``
     - FASTA format, gzipped
   * - Multi-species Alignment
     - ``resources/alignment/``
     - MAF format (EPO alignments)
   * - Population VCF
     - ``resources/pop-level-VCF/``
     - VCF with population variants (AF ≥ 0.9)
   * - Phylogenetic Tree
     - ``resources/``
     - Newick format (.nwk)

Optional Inputs
~~~~~~~~~~~~~~~

* **Gene Annotation** (GFF3/GTF) - for visualization in summary report

See documentation for detailed data preparation instructions.

Output
------

Main Output
~~~~~~~~~~~

**PHRED-scaled CADD scores**: ``results/cadd_scores/chr*.tsv.gz``

* Tab-separated, bgzip-compressed, tabix-indexed
* Format: ``CHROM  POS  REF  ALT  CADD_SCORE``
* Higher PHRED scores = More deleterious

Intermediate Outputs
~~~~~~~~~~~~~~~~~~~~

* ``results/ancestral_seq/`` - Reconstructed ancestor sequences
* ``results/derived_variants/`` - Real variants from population (benign training set)
* ``results/simulated_variants/`` - Simulated variants (deleterious-enriched training set)
* ``results/annotations/`` - VEP, PHAST, GERP annotations
* ``results/models/`` - Trained logistic regression classifier
* ``output/*.html`` - Summary reports

Future Development
------------------

High Priority
~~~~~~~~~~~~~

* Complete alignment-free branch for poorly-annotated species
* Variant derivation without population VCF using alternative methods
* Merge standard and alignment-free pipelines into unified workflow
* Add validation set support for model evaluation

Pipeline Improvements
~~~~~~~~~~~~~~~~~~~~~

* Clean up and refactor scripts for better maintainability
* Update gene annotation handling for visualization
* Production-ready cluster profiles
* Checkpoint/resume functionality for long-running jobs
* Dynamic resource allocation based on data size

Feature Enhancements
~~~~~~~~~~~~~~~~~~~~

* Enhanced indel scoring
* Integration with ab initio gene prediction (Augustus)
* Multi-species preset configurations

In a galaxy far far away, maybe some day: 
* Web interface for score queries
* Downstream analysis tools (variant prioritization, cohort analysis)

Citation
--------

If you use omniCADD, please cite:

**Original CADD Framework**:

Kircher, M., Witten, D., Jain, P. et al. A general framework for estimating the relative pathogenicity of human genetic variants. 
*Nat Genet* **46**, 310–315 (2014). https://doi.org/10.1038/ng.2892

**Species-Specific Extensions**:

* Groß, C., de Ridder, D. & Reinders, M. Predicting variant deleteriousness in non-human species: applying the CADD approach in mouse. 
  *BMC Bioinformatics* **19**, 373 (2018). https://doi.org/10.1186/s12859-018-2337-5

* Groß, C., Derks, M., Megens, HJ. et al. pCADD: SNV prioritisation in Sus scrofa. 
  *Genet Sel Evol* **52**, 4 (2020). https://doi.org/10.1186/s12711-020-0528-9

* Groß C, Bortoluzzi C, de Ridder D, Megens HJ, Groenen MAM, et al. Prioritizing sequence variants in conserved non-coding elements in the chicken genome using chCADD. 
  *PLOS Genetics* 16(9): e1009027 (2020). https://doi.org/10.1371/journal.pgen.1009027

**Ensembl VEP**:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. 
The Ensembl Variant Effect Predictor. *Genome Biology* 17(1):122 (2016). 
https://doi.org/10.1186/s13059-016-0974-4
