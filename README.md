## Repository for "A dataset of transcriptomic effects of camptothecin treatment on early zebrafish embryos"

A repository with the code, data and figures for the manuscript under consideration at "Data in Brief".
This repository contains the following parts 

1. RNA-seq mapping protocol description.
2. R code for the differential expression quantification by DESeq2 as well as for the full dataset heatmap and MA plot.
3. R code for plotting barplots of Gene Ontology Biological Process and KEGG pathway terms.
4. R code for identifying and Interferon-regulated genes and p53 target genes followed by plotting their expression heatmaps.


### 1. RNA-seq mapping protocol description.

For this study we performed mapping of RNA-sequencing reads using a new transcriptome annotation developed by the lab of Dr. Nathan Lawson.
Please see the following link [Lawson Zebrafish Transcriptome](https://www.umassmed.edu/lawson-lab/reagents/zebrafish-transcriptome/) or their [paper](https://elifesciences.org/articles/55792) for further details.

The first step is to build an index for the STAR aligner using the command in the block below. Alternatively, one can download a ready-made index the Lawson lab page mentioned above.

```
# 1. Build STAR index

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_index_v43 --genomeFastaFiles release106/Danio_rerio.GRCz11.dna.primary_assembly.fa --sjdbGTFfile v4.3.2.gtf --sjdbOverhang 149

```

The next step is to perform a first-pass alignment of reads to the genome with STAR:
```
#2. Do first-pass alignment to the genome with STAR. An example is shown. Equivalent lines need to be run in bash for all samples, ideally as a script.

STAR --runThreadN 8 --genomeDir star_index_v43  --outFileNamePrefix wt-1_  --readFilesIn 00_fastq/wt-1_R1_001.fastq 00_fastq/wt-1_R2_001.fastq
```

We next performed a second-pass alignment using some of the data from the first-pass alignment:
```
3. 2nd pass STAR alignments
# $1 == path to directory where STAR index lives
# $2 == space separated list of all splice site *tab files generated from 1-st pass
# $3 == a prefix for the STAR output, typically includes sample name
# $4 == R1 fastq file
# $5 == R2 fastq file

STAR --runThreadN 8 --genomeDir $1 --quantMode TranscriptomeSAM --sjdbFileChrStartEnd $2  --outFileNamePrefix $3 --readFilesIn $4 $5

STAR --runThreadN 8 --genomeDir star_index_v43 --quantMode TranscriptomeSAM  --sjdbFileChrStartEnd wt-1_SJ.out.tab --outFileNamePrefix wt-1_2nd_pass_ --readFilesIn 00_fastq/wt-1_R1_001.fastq 00_fastq/wt-1_R2_001.fastq
```

Before gene expression can be quantified, we need to generate an index for RSEM software. This step needs to be done once:
```
# 4. Generate an RSEM index.
rsem-prepare-reference --gtf v4.3.2.gtf release106/Danio_rerio.GRCz11.dna.primary_assembly.fa rsem_index/rsem_v43
```

The next step is to run RSEM-based gene expression quantification such as for this example file:

```
# 5. Quantify gene expression by RSEM.
rsem-calculate-expression --num-threads 8 --paired-end --alignments wt-1_2nd_passAligned.toTranscriptome.out.bam rsem_index/rsem_v43 wt-1_2nd_pass_
```

Finally, the data from individual samples need to be organized into complete data frames. This can be done using [process_data script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/R_utils/process_data.R).


### 2. R code for the differential expression quantification by DESeq2 as well as for the full dataset heatmap and MA plot.
The code for this step is included in the [deseq2_workflow script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/1_DiffExpr_code/deseq2_workflow.R). All the necessary input files for the script as well as most output files are stored in the same folder. This script also generates plots for Figure 1.

### 3. R code for plotting barplots of Gene Ontology Biological Process and KEGG pathway terms.
The code for this step is included in the [GOBP_plots script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/2_GOBP-KEGG_plots/GOBP_plots.R) and [KEGG_plots script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/2_GOBP-KEGG_plots/KEGG_plots.R). All the necessary input files for the script as well as most output files are stored in the same folder. These scripts generate plots for Figure 2.

### 4. R code for identifying and Interferon-regulated genes and p53 target genes followed by plotting their expression heatmaps.
The code for this step is included in the [IFN_p53_targets_script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/3_IRG-p53-targets_plots/IFN_p53_targets_script.R). To map known human p53 targets genes, we used [mapping_human_p53_genes_to_fish script](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper/tree/main/3_IRG-p53-targets_plots/mapping_human_p53_genes_to_fish.R)  All the necessary input files for the script as well as most output files are stored in the same folder. The plots generated by this script can be exported to pdf files used for Figures 3 and 4.







