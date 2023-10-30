# README
This is the collection of scripts that derives variants in the reference species based on the extracted ancestral sequence. First, individual level vcf files (provided by the user) is used to calculate population frequencies, where the alternative allele has a frequency above 90%. Then this information is used in the deriving of variants. Then, the simulated variants are downsampled (randomly) to match the number of derived variants.

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

Dependencies:
- bcftools
- (conda)

Python-specific dependencies are all exported in the conda environment `cadd.yml`. The pipeline can be run within this environment.

1. `generate_frequencies.py`
  Usage: `python <script.py> -v <vcf.gz file> -c <chromosome list>`
  This script take an individual level vcf file and extact the variants, together with their frequencies, if they are above 90% (0.90)

2. `derive_variants.py`
  Usage: `python derive_variants.py.py -c <chr num> -a <path to ancestor seq> -g <path to genome> -f <path to frequency files> -s <start position of region>`
  This script uses the frequency file, the reference sequence and and the extracted ancestor sequence, to extract derived variants in the reference genome, based on the criteria above.

3. `trim_variants`
  Usage: `python <script.py> -s <path to simulated variants> -d <path to derived variants> -p <prefix of simulated variant file(s)> -q <prefix of derived variant file(s)>`
  this script trims the simulated variants (randomly downsamples) down to the same amout as the derived variants.


These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
`snakemake -c4 --snakefile Snakemake_derive.sn`
