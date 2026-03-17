Output Data
===========

Complete guide to understanding and using omniCADD outputs.

Output Directory Structure
---------------------------

After pipeline completion, you'll have:

.. code-block:: text

   results/
   ├── ancestral_seq/          # Reconstructed ancestral sequences
   ├── derived_variants/       # Real variants from population (benign)
   ├── simulated_variants/     # Simulated variants (deleterious-enriched)
   ├── annotations/            # VEP, PHAST, GERP annotations
   ├── dataset/                # Combined feature matrices
   ├── models/                 # Trained classifier
   ├── cadd_scores/            # Final CADD scores ⭐
   └── figures/                # QC plots and visualizations
   
   output/
   └── *.html                  # Summary reports

Primary Output: CADD Scores
----------------------------

Location
~~~~~~~~

``results/cadd_scores/chr*.tsv.gz``

Format
~~~~~~

Tab-separated values, bgzip-compressed, tabix-indexed:

.. code-block:: text

   CHROM  POS      REF  ALT  CADD_SCORE
   1      1000000  A    G    12.5
   1      1000100  C    T    25.3
   1      1000200  G    A    8.2
   1      1000300  T    C    15.7

**Columns**:
* ``CHROM``: Chromosome
* ``POS``: Position (1-based)
* ``REF``: Reference allele
* ``ALT``: Alternate allele
* ``CADD_SCORE``: PHRED-scaled deleteriousness score

PHRED Score Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

CADD scores are PHRED-scaled, similar to quality scores:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50
   
   * - PHRED Score
     - Percentile
     - Interpretation
   * - 10
     - Top 10%
     - Moderately deleterious
   * - 15
     - Top 3%
     - Likely deleterious
   * - 20
     - Top 1%
     - Deleterious
   * - 25
     - Top 0.3%
     - Highly deleterious
   * - 30
     - Top 0.1%
     - Very highly deleterious
   * - 40+
     - Top 0.01%
     - Among most deleterious

**PHRED Calculation**:

.. code-block:: text

   PHRED = -10 * log10(1 - rank/total)

Where rank is the variant's position when sorted by raw deleteriousness score.

Using CADD Scores
-----------------

Query Specific Variants
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Single position
   tabix results/cadd_scores/chr1.tsv.gz 1:1000000-1000000
   
   # Region
   tabix results/cadd_scores/chr1.tsv.gz 1:1000000-2000000
   
   # Multiple regions from file
   tabix -R regions.bed results/cadd_scores/chr1.tsv.gz

Filter by Score Threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Extract high-impact variants (CADD > 20)
   zcat results/cadd_scores/chr1.tsv.gz | awk '$5 > 20' > high_impact.tsv
   
   # Count variants by threshold
   for threshold in 10 15 20 25 30; do
       count=$(zcat results/cadd_scores/chr*.tsv.gz | awk -v t=$threshold '$5 > t' | wc -l)
       echo "CADD > $threshold: $count variants"
   done

Annotate Your VCF
~~~~~~~~~~~~~~~~~

Add CADD scores to your variant callset:

.. code-block:: bash

   # Create VCF header
   echo '##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD deleteriousness score">' > header.txt
   
   # Annotate VCF
   bcftools annotate -a results/cadd_scores/chr1.tsv.gz \
     -c CHROM,POS,REF,ALT,CADD \
     -h header.txt \
     input.vcf.gz -Oz -o annotated.vcf.gz

Score Distribution
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Generate score histogram
   zcat results/cadd_scores/chr*.tsv.gz | \
     awk '{print $5}' | \
     sort -n | \
     uniq -c | \
     awk '{print $2"\t"$1}' > score_distribution.txt

Intermediate Outputs
--------------------

Ancestral Sequences
~~~~~~~~~~~~~~~~~~~

**Location**: ``results/ancestral_seq/{ancestor_name}/``

**Format**: FASTA, one file per chromosome

.. code-block:: text

   >1
   ATCGATCGATCGATCG...

**Use**: Reference for deriving and simulating variants

Derived Variants
~~~~~~~~~~~~~~~~

**Location**: ``results/derived_variants/singletons/``

**Format**: VCF, one file per chromosome

**Columns**: Standard VCF format

**Use**: Negative training examples (benign class)

**Count**: Typically 100k-5M variants depending on genome and population

Simulated Variants
~~~~~~~~~~~~~~~~~~

**Location**: ``results/simulated_variants/trimmed_snps/``

**Format**: VCF, one file per chromosome

**Use**: Positive training examples (deleterious class)

**Count**: Typically 500k-20M variants (5x overestimation)

Annotations
~~~~~~~~~~~

**VEP Annotations**:

.. code-block:: text

   results/annotations/vep/derived/chr*.vcf.gz
   results/annotations/vep/simulated/chr*.vcf.gz

**PHAST Conservation**:

.. code-block:: text

   results/annotations/phast/phastCons/chr*.bed
   results/annotations/phast/phyloP/chr*.bed

**GERP Scores**:

.. code-block:: text

   results/annotations/gerp/chr*.tsv

Feature Matrices
~~~~~~~~~~~~~~~~

**Location**: ``results/dataset/``

**Format**: Numpy NPZ (sparse matrices)

.. code-block:: text

   results/dataset/derived/chr1.npz
   results/dataset/simulated/chr1.npz
   results/dataset/imputation_dict.txt

**Contents**:
* Sparse feature matrices (samples × features)
* Feature names
* Sample labels

**Load with Python**:

.. code-block:: python

   import numpy as np
   from scipy.sparse import load_npz
   
   # Load feature matrix
   data = np.load('results/dataset/derived/chr1.npz', allow_pickle=True)
   features = load_npz('results/dataset/derived/chr1_features.npz')
   labels = data['labels']

Model Output
------------

Trained Model
~~~~~~~~~~~~~

**Location**: ``results/models/``

**Files**:
* ``final_model.pkl`` - Trained logistic regression classifier
* ``feature_scaler.pkl`` - Feature standardization parameters
* ``feature_coefficients.tsv`` - Feature weights

**Load Model**:

.. code-block:: python

   import joblib
   
   model = joblib.load('results/models/final_model.pkl')
   scaler = joblib.load('results/models/feature_scaler.pkl')
   
   # Use model
   probabilities = model.predict_proba(X_test)

Model Evaluation
~~~~~~~~~~~~~~~~

**Location**: ``results/models/evaluation_summary.txt``

**Contents**:

.. code-block:: text

   Cross-Validation Results:
   ========================
   Fold 1: AUC=0.89, Accuracy=0.85
   Fold 2: AUC=0.91, Accuracy=0.87
   Fold 3: AUC=0.88, Accuracy=0.84
   Fold 4: AUC=0.90, Accuracy=0.86
   Fold 5: AUC=0.89, Accuracy=0.85
   
   Mean AUC: 0.894 ± 0.011
   
   Best Parameters:
   ================
   C: 1.0
   max_iter: 2000

Feature Coefficients
~~~~~~~~~~~~~~~~~~~~

**Location**: ``results/models/feature_coefficients.tsv``

**Format**:

.. code-block:: text

   Feature                    Coefficient
   VEP_consequence_missense    0.45
   PHAST_conservation          0.38
   GERP_score                  0.32
   VEP_SIFT_score             -0.28
   ...

**Interpretation**:
* Positive coefficient = associated with deleteriousness
* Negative coefficient = associated with benign
* Magnitude = feature importance

Summary Reports
---------------

HTML Visualization
~~~~~~~~~~~~~~~~~~

**Location**: ``output/{species_comparison}.html``

**Contents**:
* Variant statistics (counts, Ti/Tv ratios)
* Interactive genome ideogram
* Phylogenetic tree visualization
* CDS coverage (if annotation provided)
* Summary tables

**Open in Browser**:

.. code-block:: bash

   # macOS
   open output/pig_cow.html
   
   # Linux
   firefox output/pig_cow.html
   
   # Or transfer to local machine
   scp cluster:path/to/output/pig_cow.html .

Quality Control Figures
~~~~~~~~~~~~~~~~~~~~~~~~

**Location**: ``results/figures/``

**Plots**:
* ``column_analysis/relevance.tsv`` - Feature importance
* ``score_distributions/`` - CADD score distributions
* ``training_curves/`` - Model training progress

Export and Share
----------------

Compress Final Scores
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create distributable archive
   tar -czf cadd_scores_pig_v1.tar.gz results/cadd_scores/
   
   # With MD5 checksum
   md5sum cadd_scores_pig_v1.tar.gz > cadd_scores_pig_v1.tar.gz.md5

Output File Sizes
-----------------

Typical file sizes for a mammalian genome:

.. list-table::
   :header-rows: 1
   
   * - Output
     - Per Chromosome
     - Whole Genome
   * - Ancestral sequence
     - 50-150 MB
     - 1-3 GB
   * - Derived variants
     - 5-20 MB
     - 100-500 MB
   * - Simulated variants
     - 20-100 MB
     - 500 MB-2 GB
   * - VEP annotations
     - 10-50 MB
     - 200 MB-1 GB
   * - PHAST scores
     - 100-300 MB
     - 2-6 GB
   * - GERP scores
     - 100-300 MB
     - 2-6 GB
   * - Feature matrices
     - 200-500 MB
     - 4-10 GB
   * - **CADD scores**
     - **500 MB-2 GB**
     - **10-40 GB**
   * - Model files
     - N/A
     - 10-100 MB

**Total storage**: 50-150 GB for complete pipeline outputs

Data Retention
--------------

Essential Files (Keep)
~~~~~~~~~~~~~~~~~~~~~~

* ✅ ``results/cadd_scores/`` - Final scores
* ✅ ``results/models/`` - Trained model
* ✅ ``output/*.html`` - Summary reports
* ✅ ``results/dataset/imputation_dict.txt`` - For reproducibility

Can Be Deleted
~~~~~~~~~~~~~~

* ⚪ ``results/ancestral_seq/`` - Can regenerate from alignment
* ⚪ ``results/derived_variants/`` - Can regenerate from VCF
* ⚪ ``results/simulated_variants/`` - Can regenerate
* ⚪ ``results/annotations/`` - Can regenerate (slow)
* ⚪ ``results/dataset/`` - Can regenerate from annotations

Cleanup Command
~~~~~~~~~~~~~~~

.. code-block:: bash

   # Remove intermediate files (saves ~70% storage)
   rm -rf results/ancestral_seq/
   rm -rf results/derived_variants/
   rm -rf results/simulated_variants/
   rm -rf results/annotations/
   
   # Keep only final outputs
   # results/cadd_scores/ (~20 GB)
   # results/models/ (~100 MB)
   # output/ (~10 MB)

Citing Outputs
--------------

When publishing results using omniCADD scores:

**Method Description**:

   "Variant deleteriousness was assessed using omniCADD scores [cite this repo], 
   a machine learning framework based on the CADD methodology (Kircher et al. 2014) 
   adapted for [species]. CADD scores were computed using functional annotations 
   from Ensembl VEP (McLaren et al. 2016) and conservation metrics from PHAST and GERP."

**Data Availability**:

   "Genome-wide CADD scores for [species] are available at [DOI/URL]."

See Also
--------

* :doc:`workflow` - How outputs are generated
