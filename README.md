# dmvar
Identify DNA sequence variants in Drosophila mutant lines

The key thing this pipeline introduces is the identification of variants which are novel in a particular mutant fly line compared to the parent and other progeny of that parent (the cohort). Depending on the line the variants are filtered to those which are heterozygous (when balanced due to lethality) or homozygous, as specified in the samplesheet. The location of variants of interest is specified by chromosome or chromosome arm, again in the samplesheet.

Input - REQUIRED: nf-core/sarek samplesheet, dmvar VCF samplesheet, fastq files

OPTIONAL: Gene summary file (e.g. http://ftp.flybase.org/releases/FB2021_06/precomputed_files/genes/automated_gene_summaries.tsv.gz), gene map file (e.g. http://ftp.flybase.org/releases/FB2021_06/precomputed_files/genes/gene_map_table_fb_2021_06.tsv.gz)

Output - Excel files of filtered variants per cohort, per mutant (plus bam files of mapped reads, combined, genotyped GVCF and VCFs of filtered SNPs and indels)

Variants are considered only if they are homozygous in the parent and either homozygous (but non-parental) in a single progeny line or heterozygous alt in a single progeny line

## Running dmvar

1. Run nf-core/sarek to get GVCFs per sample 

`nextflow run nf-core/sarek -r 2.7.1 --input <samplesheet> -c gurdon.config --genome dm6 --save_bam_mapped --tools haplotypecaller --generate_gvcf --outdir <outdir>`

2. Run dmvar nextflow script to combine GVCFs, call variants, filter and annotate them

`nextflow run dmvar.nf -c gurdon.config -with-singularity nfcore-sarek-2.7.1.img --ss <dmvar_samplesheet> --ref dm6.fa --dict dm6.dict --fai dm6.fa.fai --snpeff_ref BDGP6.28.99 --outdir dmvar_outdir`

3. Run python script to generate variants of interest in Excel files per parental control

N.b. This requires bcftools and python imports subprocess, pandas and openpyxl

`python ~/code_development/dmvar/dmvar_results.py -s vcf_samplesheet.csv -v combined_filtered_pass.vcf.snpeff.vcf.gz -f HIGH,MODERATE -g ~/code_development/dmvar/automated_gene_summaries.tsv -m ~/code_development/dmvar/gene_map_table_fb_2021_06.tsv.gz`

## Samplesheets

Example samplesheet for sarek (no header required): `sarek_samplesheet_example.tsv`

Example samplesheet VCF for dmvar (header required): `dmvar_samplesheet_example.csv`

- Parents should have 'None' in *chromosome*, *type* and *control* columns
- *chromosome* can be a stem e.g. chr2 will result in variants from chr2L and chr2R being included (this makes sense for Drosophila, but might not for other genomes)

## Accessory files

gene_map_table_fb_2021_06.tsv.gz - File downloaded from Flybase describing the locations of genes. Used to annotate the results

automated_gene_summaries.tsv - File downloaded from Flybase describing known functions etc. of genes. Used to annotate the results


## Requires nf-core singularity images

depot.galaxyproject.org-singularity-snpeff-5.0--hdfd78af_1.img (nfcore-sarek-2.7.1.img contains an earlier version of snpEff which is not compatible)

nfcore-sarek-2.7.1.img

## Config files

nextflow.config - indicates that Nextflow should use the snpeff Singularity image to run the annotate_variants process (will need to be edited for other systems)

gurdon.config - describes how nextflow should run through Slurm on the Gurdon cluster 
