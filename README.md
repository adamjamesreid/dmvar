# dmvar
Identify DNA sequence variants in Drosophila mutant lines

## Requires nf-core singularity images

depot.galaxyproject.org-singularity-snpeff-5.0--hdfd78af_1.img

nfcore-sarek-2.7.1.img

## Running dmvar

1. Run nf-core/sarek to get GVCFs per sample 

`nextflow run nf-core/sarek -r 2.7.1 --input <samplesheet> -c gurdon.config --genome dm6 --save_bam_mapped --tools haplotypecaller --generate_gvcf --outdir <outdir>`

2. Run dmvar nextflow script to combine GVCFs, call variants, filter and annotate them

`nextflow run dmvar.nf -c gurdon.config -c nextflow.config -with-singularity nfcore-sarek-2.7.1.img --ss <dmvar_samplesheet> --ref dm6.fa --dict dm6.dict --fai dm6.fa.fai --snpeff_ref BDGP6.28.99 --outdir dmvar_outdir`

3. Run python script to identify variants of interest in Excel files per parental control

N.b. This incorporates remove_parental2.py and currently has a hidden requirement for bcftools and python imports subprocess, pandas and openpyxl! Genesummary data is also not optional, but should be made so.

`python ~/code_development/dmvar/dmvar_results.py -g ~/code_development/dmvar/automated_gene_summaries.tsv ../NVS024_vcf_samplesheet.csv -r ~/code_development/dmvar/remove_parental2.py -v snps_filtered_pass.vcf.snpeff.vcf`

## Samplesheets

Example samplesheet for sarek (no header required): sarek_samplesheet_example.tsv

Example samplesheet VCF for dmvar (header required): dmvar_samplesheet_example.csv

Parents should have 'None' in chromosome, type and control columns
