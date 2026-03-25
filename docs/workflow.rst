Workflow Overview
=================

omniCADD processes genomic data through 8 sequential steps to generate CADD scores.

Pipeline Architecture
---------------------

The pipeline follows this general flow:

.. code-block:: text

   External Resources (Genome, Alignment, VCF, Tree)
                    ↓
         Step 1: Extract Ancestral Sequences
                    ↓
      ┌─────────────┴─────────────┐
      ↓                           ↓
   Step 2: Derive           Step 3: Simulate
   Variants (Benign)        Variants (Deleterious)
      ↓                           ↓
      └─────────────┬─────────────┘
                    ↓
         Step 5: Annotate Variants
         (VEP, PHAST, GERP)
                    ↓
         Step 6: Combine Annotations
         (Feature Matrices)
                    ↓
         Step 7: Train Model
         (Logistic Regression)
                    ↓
         Step 8: Score All Variants
                    ↓
              CADD Scores

   Step 4: Summary Report (runs independently after Steps 1-3)

Workflow Steps
--------------

Step 1: Extract Ancestral Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Reconstruct the most recent common ancestor (MRCA) sequence from multi-species alignments.

**Input**:
  * Multi-species alignment (MAF format)
  * Phylogenetic tree (Newick format)
  * Reference genome

**Process**:
  1. Clean MAF alignment files (remove ambiguous nucleotides)
  2. Identify ancestral node between species of interest and comparison species
  3. Deduplicate sequences and reorder species
  4. Correct strand orientation
  5. Sort alignment blocks by chromosome position
  6. Extract ancestor sequences in reference coordinates

**Output**:
  * ``results/ancestral_seq/{ancestor_name}/chr*.fa``

**Key Scripts**:
  * ``workflow/scripts/step_1_extract_ancestor/*.py``

**Runtime**: ~2-6 hours for mammalian genome

---

Step 2: Derive Variants
~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Identify real variants from population VCF that differentiate the reference from the ancestor.

**Input**:
  * Population VCF with allele frequencies
  * Reference genome
  * Ancestral sequences (from Step 1)

**Process**:
  1. Extract high-frequency variants (AF ≥ 0.9) from population
  2. Apply 5-case logic comparing ancestor, reference, and population alleles:
     
     * Case 1: Ref = Anc, Pop has derived variant
     * Case 2: Ref ≠ Anc, Ref is fixed in population
     * Case 3: Multi-allelic with derived fixation
     * Case 4: Discordant Ref vs Pop
     * Case 5: Complex sites
     
  3. Filter for singleton SNPs (remove clustered variants)

**Output**:
  * ``results/derived_variants/singletons/chr*.vcf``

**Key Scripts**:
  * ``workflow/scripts/step_2_derive_variants/*.py``

**Runtime**: ~1-3 hours

**Note**: Requires population VCF with AF ≥ 0.9 variants

---

Step 3: Simulate Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Generate simulated variants enriched for deleterious mutations for training.

**Input**:
  * Reference genome
  * Ancestral sequences
  * Mutation parameter estimates

**Process**:
  1. Compute substitution rates from ancestor ↔ reference differences
  2. Scale parameters by overestimation factor (default 5x)
  3. Simulate SNPs per chromosome using mutation model
  4. Optionally simulate INDELs (computationally intensive)
  5. Filter variants on gaps or missing ancestral sequence

**Output**:
  * ``results/simulated_variants/trimmed_{snps|indels}/chr*.vcf``

**Key Scripts**:
  * ``workflow/scripts/step_3_simulate_variants/*.py``

**Algorithm**: Based on Kircher et al. 2014 original implementation

**Runtime**: ~2-4 hours for SNPs, +6-12 hours for INDELs

---

Step 4: Summary Report
~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Generate interactive HTML visualization of genome and variant properties.

**Input**:
  * Derived and simulated variants
  * Ancestral sequences
  * Phylogenetic tree
  * Optional: Gene annotation (GFF3/GTF)

**Output**:
  * ``output/{species_comparison}.html``

**Contents**:
  * Variant statistics (counts, Ti/Tv ratios)
  * Ideogram-style chromosome browser
  * Phylogenetic tree visualization
  * CDS coverage plots (if annotation provided)
  * Interactive tables and plots

**Key Scripts**:
  * ``workflow/scripts/step_4_summary_report/*.Rmd``

**Technology**: R with rmarkdown, ggplot2, ggtree

**Runtime**: ~30 minutes - 2 hours

**Note**: Can run in Docker container if R packages unavailable

---

Step 5: Annotate Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Apply functional and conservation annotations to all variants.

**Input**:
  * Derived variants (from Step 2)
  * Simulated variants (from Step 3)
  * Reference genome
  * Multi-species alignments

**Annotation Types**:

**A. Functional Annotation (VEP/SNPEff)**:
  * Variant consequence prediction (missense, nonsense, synonymous, etc.)
  * SIFT deleteriousness scores (when available)
  * Grantham matrix for amino acid substitution severity
  * Gene context and transcript information

**B. Conservation Annotation (PHAST)**:
  * ``phyloFit``: Train phylogenetic substitution models
  * ``phastCons``: Conservation scores (0-1 scale)
  * ``phyloP``: P-values for accelerated/conserved evolution

**C. Conservation Annotation (GERP)**:
  * ``gerpcol``: Compute expected vs observed substitution rates
  * GERP rejection scores (higher = more conserved)
  * Element-based conservation scores

**Output**:
  * ``results/annotations/vep/{derived|simulated}/chr*.vcf.gz``
  * ``results/annotations/phast/*.bed``
  * ``results/annotations/gerp/*.tsv``

**Key Scripts**:
  * ``workflow/scripts/step_5_annotate_variants/*.py``

**Runtime**: 
  * VEP: ~4-8 hours
  * PHAST: ~6-12 hours (parallel across alignment chunks)
  * GERP: ~8-16 hours

---

Step 6: Combine Annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Merge all annotation types into unified feature matrices for model training.

**Input**:
  * VEP annotations (from Step 5)
  * PHAST scores (from Step 5)
  * GERP scores (from Step 5)
  * Multi-species alignments

**Process**:
  1. Split alignments into manageable chunks (~10k positions)
  2. Convert MAF → FASTA format for conservation tools
  3. Prune alignment columns with excessive gaps (>50%)
  4. Run GERP and PHAST on alignment chunks in parallel
  5. Merge VEP + GERP + PHAST outputs by genomic position
  6. Create sparse numpy NPZ datasets with feature matrices
  7. Calculate imputation means for missing values
  8. Analyze feature relevance and distributions

**Output**:
  * ``results/dataset/{derived|simulated}/chr*.npz`` (sparse matrices)
  * ``results/dataset/imputation_dict.txt`` (feature means for missing values)
  * ``results/figures/column_analysis/relevance.tsv`` (feature importance)

**Key Scripts**:
  * ``workflow/scripts/step_6_combine_annotations/*.py``

**Runtime**: ~8-16 hours (highly parallel)

**Features Generated**: 50-100+ features depending on annotation configuration

---

Step 7: Train & Test Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Build machine learning classifier to predict variant deleteriousness.

**Input**:
  * Combined feature matrices (from Step 6)
  * Labels: Derived = 0 (benign), Simulated = 1 (deleterious)

**Process**:
  1. **Fold Data**: Split combined dataset into K folds (default 5)
  2. **Grid Search**: K-fold cross-validation across hyperparameter space
     
     * Regularization C: [0.001, 0.01, 0.1, 1, 10]
     * Max iterations: [1000, 2000, 5000]
     
  3. **Train Models**: Train logistic regression for each parameter combination
  4. **Evaluate**: Calculate AUC, accuracy, precision, recall per fold
  5. **Final Model**: Train on all data with best hyperparameters
  6. **Export**: Save model, scaler, and feature coefficients

**Output**:
  * ``results/models/final_model.pkl`` (trained classifier)
  * ``results/models/feature_scaler.pkl`` (standardization parameters)
  * ``results/models/feature_coefficients.tsv`` (feature weights)
  * ``results/models/evaluation_summary.txt`` (cross-validation metrics)

**Algorithm**: scikit-learn LogisticRegression with L2 regularization

**Key Scripts**:
  * ``workflow/scripts/step_7_train_test_model/*.py``

**Runtime**: ~2-6 hours (depends on dataset size and grid search)

**Typical Performance**: AUC 0.85-0.95 on held-out test sets

---

Step 8: Score Variants
~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Generate genome-wide PHRED-scaled CADD scores for all possible variants.

**Input**:
  * Trained model and scaler (from Step 7)
  * Reference genome
  * Conservation annotations (PHAST, GERP)

**Process**:
  1. Generate all possible SNPs per chromosome in blocks (~10k positions/block)
  2. Annotate all variants with VEP
  3. Intersect with conservation scores (PHAST, GERP)
  4. Prepare feature matrices matching training format
  5. Apply trained model to predict deleteriousness probabilities
  6. Rank variants and convert to PHRED scale:
     
     ``PHRED = -10 * log10(1 - rank/total)``
     
  7. Merge blocks, sort by chromosome:position
  8. Compress with bgzip and index with tabix

**Output**:
  * ``results/cadd_scores/chr*.tsv.gz`` (bgzip-compressed, tabix-indexed)
  * ``results/cadd_scores/scoring_summary.txt`` (summary statistics)

**Format**:
  .. code-block:: text
  
     CHROM  POS      REF  ALT  CADD_SCORE
     1      1000000  A    G    12.5
     1      1000100  C    T    25.3

**Key Scripts**:
  * ``workflow/scripts/step_8_score_variants/*.py``

**Runtime**: ~12-48 hours (entire genome, all possible SNPs)

**PHRED Score Interpretation**:
  * PHRED 10: Top 10% most deleterious
  * PHRED 20: Top 1% most deleterious
  * PHRED 30: Top 0.1% most deleterious

---

Parallelization Strategy
-------------------------

The pipeline is designed for efficient parallel execution:

Chromosome-Level Parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most steps process chromosomes independently:

.. code-block:: bash

   # Snakemake automatically parallelizes across chromosomes
   snakemake --cores 16
   
   # 16 cores can process ~8-16 chromosomes in parallel

Within-Step Parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some steps have internal parallelization:

* **Step 6**: Alignment chunks processed in parallel
* **Step 7**: Cross-validation folds run in parallel
* **Step 8**: Variant blocks scored in parallel

Resource Requirements by Step
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 15 35
   
   * - Step
     - CPU
     - RAM
     - Runtime
     - Bottleneck
   * - 1. Extract Ancestor
     - 2-4
     - 8 GB
     - 2-6 h
     - I/O (large MAF files)
   * - 2. Derive Variants
     - 1-2
     - 4 GB
     - 1-3 h
     - VCF parsing
   * - 3. Simulate Variants
     - 1-2
     - 4 GB
     - 2-4 h
     - Random number generation
   * - 4. Summary Report
     - 2-4
     - 8 GB
     - 0.5-2 h
     - R package loading
   * - 5. Annotate (VEP)
     - 4-8
     - 16 GB
     - 4-8 h
     - Database lookups
   * - 5. Annotate (PHAST)
     - 8-16
     - 32 GB
     - 6-12 h
     - Phylogenetic models
   * - 5. Annotate (GERP)
     - 8-16
     - 32 GB
     - 8-16 h
     - Substitution rate calc
   * - 6. Combine
     - 8-16
     - 64 GB
     - 8-16 h
     - Sparse matrix ops
   * - 7. Train Model
     - 8-16
     - 32 GB
     - 2-6 h
     - Grid search
   * - 8. Score Variants
     - 16-32
     - 64 GB
     - 12-48 h
     - Prediction on all SNPs

Checkpointing and Resume
-------------------------

Snakemake automatically handles checkpointing:

.. code-block:: bash

   # If pipeline is interrupted
   snakemake --use-conda --cores 16
   
   # Snakemake resumes from last completed step

Force Re-run Specific Steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Re-run from Step 6 onwards
   snakemake --use-conda --cores 16 --forcerun combine_annotations
   
   # Re-run only model training
   snakemake --use-conda --cores 16 --forcerun train_model

See Also
--------

* :doc:`configuration` - Configure pipeline parameters
