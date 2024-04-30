# README
This is the collection of scripts that derives variants in the reference species based on the extracted ancestral sequence. First, individual level vcf files (provided by the user) is used to calculate population frequencies, where the alternative allele has a frequency above 90%. Then this information is used in the deriving of variants. 

The identification is based on five criteria:
|Cases|1|2|3|4|5|
|:--|:--|:--|:--|:--|:--|
|*Ancestor* | C | C	| C	| C	| C |
|*Reference* | C | A>0.9 | A | A>0.9 | A |
|*Alternative* | A>0.9 | G | G>0.9 | C | A |

- **Case 1**: The reference genome and ancestor are the same at position n,
		but in the population data, there is a nearly fixed derived variant with >0.9
- **Case 2**: The reference genome and ancestor are different at position n,
		and the reference allele has a frequency >0.9
- **Case 3**: The reference genome and ancestor are different at position n,
		but there is a third allele in the population with a frequency >0.9
- **Case 4**: The reference genome and ancestor are different, but the alternative allele
		in the population data is not the one in the reference genome, i.e. is most
		likely ancestral. Here the population reference allele at >0.9 is what is desired
- **Case 5**: The ancestral and reference allele are different at position n

*Case X: The reference genome and ancestor are different, but there is no data on that position*
		*in the population file. If a true variant and true missingness, it can be *assumed to be 100% fixed*
		in the reference species / population, and hence be monomorphic in the reference.*
		**not yet implemented**

Python dependencies:
- biopython
- pysam
- mafTools
- vcftools
- snakemake
- (conda)

Python-specific dependencies are all exported in the conda environment `simulation.yml`. The pipeline handles user proivided `.vcf`-files. This should contain data for multiple samples in a population, from which the allele frequency should be derived. 

**To do**
- add derived variant summary stats to the stats report
- add alternative ways of deriving variants when pop. level data is not available


