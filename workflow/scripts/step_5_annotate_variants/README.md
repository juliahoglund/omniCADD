# README
This is the collection of scripts that annotated the simulate dand derived variants. First [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)<sup>[ref](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4)</sup> annotation is performed. Then the files are processed and reformatted, and variants are further annotated with [GERP](http://mendel.stanford.edu/sidowlab/downloads/gerp/index.html)<sup>[ref](https://genome.cshlp.org/content/15/7/901)</sup>, [PhastCons](http://mendel.stanford.edu/sidowlab/downloads/gerp/index.html)<sup>[ref](https://genome.cshlp.org/content/15/7/901)</sup>, and [PhyloP](http://mendel.stanford.edu/sidowlab/downloads/gerp/index.html)<sup>[ref](https://genome.cshlp.org/content/15/7/901)</sup>. Finally information about repeat position are added, and all anotations are joined and merged into one file per chromosome.

**Abbreviation used for the consequences are:**
|Abr.|Cons.|Abr.|Cons|
|:--|:--|:--|:--|
|**SG**|Stop_Gained;|**NS**|Non_Synonymous (missense);|
|**IF**|Inframe_Insertion;|**FS**|Frame_Shift;|
|**SL**|Stop_Lost;|**CS**|Canonical_Splice (splice donor);|
|**S**|Splice_Site (splice donor);|**NC**|Noncoding_Change (non-coding exon);|
|**SN**|Synonymous;|**IG**|Intergenic;|
|**DN**|Downstream;|**UP**|Upstream;|
|**R**|Regulatory_Region|**U5**|5Prime_UTR;|
|**U3**|3Prime_UTR;|**I**|Intronic;|
|**O**| Unknown|

Dependencies
- vep
- PHAST
- GERP
- SIFT
- seqtk
- biopython
- wig2bed (bedops)
- perl

Dependencies are all exported in the conda environment `annotation.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`

#### scripts adapted from other software / pipelines:

`maf2fasta.pl`: [mugsy alignment tool](https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl)
`split_alignments.py`: [the virtual laboratory](https://thevirtuallaboratory.com/blog/splitting-a-multi-fasta)
`gerp_to_position.py`: [GenErode pipeline](https://github.com/NBISweden/GenErode)
`prune_cols.py`: [compbio-utils](https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py)
`make-SIFT-db-all.pl`: [pauline-ng;SIFT4G_Create_Genomic_DB](https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB)


## TODO
- implement the unimplemented rules and annotations
2025-03-06: clean-up untested
