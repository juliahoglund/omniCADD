Input Data Preparation
======================

Complete guide to preparing all required input data for omniCADD.

Data Requirements Overview
---------------------------

omniCADD requires four types of input data:

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 35
   
   * - Data Type
     - Required?
     - Typical Size
     - Purpose
   * - Reference Genome
     - ✅ Required
     - 0.5-3 GB
     - Reference sequences
   * - Multi-species Alignment
     - ✅ Required
     - 50-500 GB
     - Ancestral reconstruction
   * - Population VCF
     - ✅ Required
     - 1-100 GB
     - Derived variants
   * - Phylogenetic Tree
     - ✅ Required
     - < 1 MB
     - Species relationships
   * - Gene Annotation
     - ⚪ Optional
     - 10-100 MB
     - Visualization only

Directory Structure
-------------------

Place all input data in the ``resources/`` directory:

.. code-block:: text

   resources/
   ├── genome/              # Reference genome FASTA
   ├── alignment/           # Multi-species alignments (MAF)
   ├── pop-level-VCF/       # Population variant data
   ├── tree_*.nwk          # Phylogenetic tree
   └── *.gff3 / *.gtf      # Optional gene annotations

1. Reference Genome
-------------------

Download Reference
~~~~~~~~~~~~~~~~~~

**Ensembl** (Recommended, when available):

.. code-block:: bash

   cd resources/genome
   
   # Pig
   wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
   
   # Chicken
   wget https://ftp.ensembl.org/pub/current_fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
   
   # Mouse
   wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz

**NCBI RefSeq**:

.. code-block:: bash

   # Navigate to NCBI FTP
   # https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
   
   # Download genomic FASTA
   wget <URL_to_genomic.fna.gz>

**UCSC**:

.. code-block:: bash

   # Example: Human
   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Validate Genome
~~~~~~~~~~~~~~~

.. code-block:: bash

   # Check file size
   ls -lh *.fa.gz
   
   # Count chromosomes
   zcat *.fa.gz | grep "^>" | wc -l
   
   # View first sequences
   zcat *.fa.gz | head -n 50

Index Genome
~~~~~~~~~~~~

The pipeline auto-indexes, but you can pre-index:

.. code-block:: bash

   samtools faidx genome.fa.gz

**Expected output**: ``genome.fa.gz.fai``

Chromosome Naming
~~~~~~~~~~~~~~~~~

Ensure consistent chromosome naming across all files:

.. list-table::
   :header-rows: 1
   
   * - Source
     - Format
     - Example
   * - Ensembl
     - ``1``, ``2``, ``X``
     - ``1``
   * - UCSC
     - ``chr1``, ``chr2``, ``chrX``
     - ``chr1``
   * - NCBI
     - ``NC_*`` or ``1``
     - ``NC_010443.5`` or ``1``

**Standardize if needed**:

.. code-block:: bash

   # Add "chr" prefix
   zcat genome.fa.gz | sed 's/^>\([0-9XYM]\)/^>chr\1/' | bgzip > genome_chr.fa.gz
   
   # Remove "chr" prefix
   zcat genome.fa.gz | sed 's/^>chr/^>/' | bgzip > genome_nochr.fa.gz

2. Multi-Species Alignment
---------------------------

Download EPO Alignment
~~~~~~~~~~~~~~~~~~~~~~

**Ensembl EPO** (Recommended):

.. code-block:: bash

   cd resources/alignment
   
   # 44-way mammals
   wget -r -nd --no-parent -e robots=off \
     -A '44_mammals.epo.*.maf.gz' \
     https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/44_mammals.epo/
   
   # 22-way primates
   wget -r -nd --no-parent -e robots=off \
     -A '22_primates.epo.*.maf.gz' \
     https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/22_primates.epo/
   
   # 12-way fish
   wget -r -nd --no-parent -e robots=off \
     -A '12_fish.epo.*.maf.gz' \
     https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/12_fish.epo/

**UCSC Multiz**:

.. code-block:: bash

   # Example: Human 100-way
   wget -r -nd --no-parent \
     https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/

Validate Alignment
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Count alignment blocks
   zcat *.maf.gz | grep "^a score" | wc -l
   
   # List species in alignment
   zcat *.maf.gz | grep "^s " | cut -f2 -d' ' | cut -f1 -d'.' | sort -u
   
   # View first block
   zcat chr1.maf.gz | head -n 50

**Expected output**:

.. code-block:: text

   a score=12345.0
   s sus_scrofa.1 100 50 + 274330532 ATCGATCG...
   s bos_taurus.5 200 50 + 158534110 ATCGATCG...
   s equus_caballus.10 300 50 + 144923456 ATCGATCG...

Alignment Quality
~~~~~~~~~~~~~~~~~

**Species Selection**:
* Minimum: 3 species (reference + 2 outgroups)
* Recommended: 10-50 species
* Optimal: Diverse phylogeny with close and distant relatives

**Quality Metrics**:
* EPO alignments: Highest quality, strand-aware
* Multiz: Good quality, widely available
* Custom: Depends on alignment tool

Storage Requirements
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   
   * - Alignment
     - Species
     - Compressed
     - Decompressed
   * - 44 mammals EPO
     - 44
     - ~280 GB
     - ~2 TB
   * - 22 primates EPO
     - 22
     - ~120 GB
     - ~800 GB
   * - 12 fish EPO
     - 12
     - ~50 GB
     - ~300 GB
   * - UCSC 100-way
     - 100
     - ~500 GB
     - ~3 TB

**Recommendation**: Keep files compressed; pipeline reads ``.maf.gz`` directly

3. Population VCF
-----------------

Prepare Population VCF
~~~~~~~~~~~~~~~~~~~~~~

If you have individual VCF files:

.. code-block:: bash

   cd resources/pop-level-VCF
   
   # Merge individual VCFs
   bcftools merge sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
     -Oz -o population_merged.vcf.gz
   
   # Calculate allele frequencies
   bcftools +fill-tags population_merged.vcf.gz \
     -Oz -o population_AF.vcf.gz -- -t AF,AN,AC
   
   # Index VCF
   tabix -p vcf population_AF.vcf.gz

Validate VCF
~~~~~~~~~~~~

.. code-block:: bash

   # Count variants
   bcftools view -H population.vcf.gz | wc -l
   
   # Check AF field
   bcftools query -f '%CHROM\t%POS\t%INFO/AF\n' population.vcf.gz | head
   
   # Verify chromosome names match genome
   bcftools view -h population.vcf.gz | grep "^##contig"

**Expected output**:

.. code-block:: text

   1       1000    0.95
   1       2000    0.87
   1       3000    1.00

Requirements
~~~~~~~~~~~~

* **Allele Frequency**: Must have AF in INFO field
* **Quality**: QUAL > 30, DP > 10 recommended
* **Format**: VCF or BCF, bgzip-compressed, tabix-indexed

Sample Size Recommendations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   
   * - Sample Size
     - Confidence
     - Use Case
   * - 10-30
     - Low
     - Pilot studies
   * - 50-100
     - Medium
     - Most species
   * - 100+
     - High
     - Ideal
   * - 1000+
     - Very High
     - Population genomics

4. Phylogenetic Tree
---------------------

Download Tree
~~~~~~~~~~~~~

.. code-block:: bash

   cd resources
   
   # Ensembl species tree (matching alignment)
   wget https://ftp.ensembl.org/pub/current_compara/species_trees/44_eutherian_mammals_EPO_default.nh
   
   # Rename appropriately
   mv 44_eutherian_mammals_EPO_default.nh tree_43_mammals.nwk

Format Requirements
~~~~~~~~~~~~~~~~~~~

**Newick format** (``.nwk``, ``.nh``, ``.tree``):

.. code-block:: text

   ((Sus_scrofa:0.01,Bos_taurus:0.02):0.05,(Equus_caballus:0.03,Canis_lupus:0.04):0.06);

**Requirements**:
* Species names must match alignment exactly (case-sensitive)
* Branch lengths in millions of years (to match gerp)
* Tree must include all species in alignment

Validate Tree
~~~~~~~~~~~~~

.. code-block:: bash

   # View tree
   cat tree_43_mammals.nwk
   
   # Validate format with Python
   python -c "from Bio import Phylo; Phylo.read('tree_43_mammals.nwk', 'newick')"
   
   # Extract species list
   grep -o "[A-Za-z_]*:" tree_43_mammals.nwk | tr -d ':' | sort

Create Custom Tree
~~~~~~~~~~~~~~~~~~

If you need a custom tree:

**Using IQTREE**:

.. code-block:: bash

   # From alignment
   iqtree -s alignment.fa -m GTR+G -bb 1000
   
   # Output: alignment.fa.treefile

**Using RAxML**:

.. code-block:: bash

   # Fast bootstrap
   raxml-ng --all --msa alignment.fa --model GTR+G --threads 8

**Using TimeTree** (divergence times):

Visit: http://www.timetree.org/
*Note*: TimeTree provides divergence times in billion years, convert to millions for gerp compatibility.

5. Gene Annotation (Optional)
------------------------------

Download Annotation
~~~~~~~~~~~~~~~~~~~

**Ensembl**:

.. code-block:: bash

   cd resources
   
   # GFF3 format
   wget https://ftp.ensembl.org/pub/current_gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gff3.gz
   
   # GTF format
   wget https://ftp.ensembl.org/pub/current_gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gtf.gz

**NCBI**:

.. code-block:: bash

   # Download from NCBI FTP
   wget <URL_to_genomic.gff.gz>

**Note**: Gene annotations are **optional** and only used for visualization in the summary report.

Data Validation Checklist
--------------------------

Before Running Pipeline
~~~~~~~~~~~~~~~~~~~~~~~

✅ **Reference Genome**:

.. code-block:: bash

   # Exists and is readable
   test -f resources/genome/*.fa.gz && echo "✓ Genome found"
   
   # Can be decompressed
   zcat resources/genome/*.fa.gz | head -n 10

✅ **Alignment**:

.. code-block:: bash

   # Files exist
   ls resources/alignment/*.maf.gz | wc -l
   
   # Contains reference species
   zcat resources/alignment/*.maf.gz | grep "^s sus_scrofa" | head -n 1

✅ **Population VCF**:

.. code-block:: bash

   # Has AF field
   bcftools view -h resources/pop-level-VCF/*.vcf.gz | grep "##INFO=<ID=AF"
   
   # Is indexed
   test -f resources/pop-level-VCF/*.vcf.gz.tbi && echo "✓ VCF indexed"

✅ **Phylogenetic Tree**:

.. code-block:: bash

   # File exists
   test -f resources/tree_*.nwk && echo "✓ Tree found"
   
   # Valid Newick format
   python -c "from Bio import Phylo; Phylo.read('resources/tree_43_mammals.nwk', 'newick')" && echo "✓ Tree valid"

✅ **Chromosome Names Match**:

.. code-block:: bash

   # Compare genome and VCF chromosome names
   zcat resources/genome/*.fa.gz | grep "^>" | head
   bcftools view -H resources/pop-level-VCF/*.vcf.gz | cut -f1 | sort -u | head
