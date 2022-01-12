#!/usr/bin/env nextflow

// input is a samplesheet defining the experiment  
// samplesheet is tricky to use though!! so I've largely just set the parameters, while getting the gvcf filenames from the sample sheet
params.ss = 'samplesheet.csv'
params.ref = "$baseDir/dm6.fa"
params.fai = "$baseDir/dm6.fa.fai"
params.dict = "$baseDir/dm6.dict"
params.outdir = "$baseDir/outdir"
params.snpeff_ref = 'BDGP6.28.99'


log.info """\
        DROSOPHILA MUTANT VARIANT CALLING PIPELINE (dmvar)
        ==================================================
        Requires cohort samplesheet with single GVCFs per sample e.g. from nf-core/sarek
        Samplesheet: 	${params.ss}
        Reference:	${params.ref}
        Ref .fai:	${params.fai}
        Ref .dict:	${params.dict}
        Output dir:	${params.outdir}
        snpEff ref:	${params.snpeff_ref}
        """
        .stripIndent()

// Check input path parameters to see if they exist
checkPathParamList = [
    params.ss, params.ref, params.fai, params.dict
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//Set up reference channels
Channel
    .fromPath(params.ref)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{ref_ch}

Channel
    .fromPath(params.fai)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{fai_ch}

Channel
    .fromPath(params.dict)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{dict_ch}

// Read in samplesheet and get vcfs
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.vcf)}
    .set{ samples_ch }

// Make a channel of vcf indexes
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.vcf + '.tbi', checkIfExists: true)}
    .set{ samples_index_ch }

process combineGVCFs {
    input:
    path ref from ref_ch
    path fai from fai_ch
    path dict from dict_ch
    path vcf from samples_ch.collect()
    path vcf_index from samples_index_ch.collect()

    output:
    path 'combined.g.vcf' into combined_ch

    script:
    """
    echo $vcf
    echo $vcf_index
    gatk CombineGVCFs -R $ref -O combined.g.vcf -V ${vcf.join(' -V ')}
    """
}

process genotypeGVCFs {
    publishDir params.outdir, mode:'copy'

    input:
    path ref from ref_ch
    path fai from fai_ch
    path dict from dict_ch
    path vcf from combined_ch

    output:
    path 'combined.GenotypeGVCFs.vcf.gz' into combined_genotyped_ch
    path 'combined.GenotypeGVCFs.vcf.gz.tbi' into combined_genotyped_index_ch

    script:
    """
    gatk GenotypeGVCFs -R $ref -V $vcf -O combined.GenotypeGVCFs.vcf
    bgzip combined.GenotypeGVCFs.vcf
    tabix -p vcf combined.GenotypeGVCFs.vcf.gz
    """
    
}

//subset SNPs and indels
process subset_variants {
    input:
    path comb_gvcf from combined_genotyped_ch
    path comb_gvcf_index from combined_genotyped_index_ch

    output:
    path 'snps.vcf.gz' into snps_ch
    path 'snps.vcf.gz.tbi' into snps_index_ch
    path 'indels.vcf.gz' into indels_ch
    path 'indels.vcf.gz.tbi' into indels_index_ch

    script:
    """
    gatk SelectVariants -V $comb_gvcf -select-type SNP -O snps.vcf.gz
    gatk SelectVariants -V $comb_gvcf -select-type INDEL -O indels.vcf.gz
    """
}

//filter SNPs
process filter_snps {
    input:
    path snp_vcf from snps_ch
    path snp_vcf_index from snps_index_ch

    output:
    path 'snps_filtered_pass.vcf' into snps_filt_ch

    script:
    """
    gatk VariantFiltration -V $snp_vcf -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR >3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O snps_filtered.vcf.gz

    bcftools view -f PASS -o snps_filtered_pass.vcf snps_filtered.vcf.gz
    """
}

//filter indels
process filter_indels {
    input:
    path indel_vcf from indels_ch
    path indel_vcf_index from indels_index_ch
    
    output:
    path 'indels_filtered_pass.vcf' into indels_filt_ch

    script:
    """
    gatk VariantFiltration -V $indel_vcf -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' -O indels_filtered.vcf.gz

    bcftools view -f PASS -o indels_filtered_pass.vcf indels_filtered.vcf.gz
    """
}

//Annotate variants with snpeff
process annotate_variants {
    publishDir params.outdir, mode:'copy'

    input:
    path vcf from snps_filt_ch.mix(indels_filt_ch) // Mix snps and indel channels together

    output:
    path "${vcf}.snpeff.vcf" into ann_vcf_ch
    //path "${vcf}.snpeff.vcf.gz.tbi" into ann_vcf_ch_index

    script:
    """
    snpEff ann ${params.snpeff_ref} -dataDir /tmp/snpEff $vcf > "${vcf}.snpeff.vcf"
    #bgzip ${vcf}.snpeff.vcf
    #tabix -o vcf ${vcf}.snpeff.vcf.gz
    """
}

//zip and index annotated VCFs
process index_vcf
{
    publishDir params.outdir, mode:'copy'

    input:
    path ann_vcf from ann_vcf_ch

    output:
    path "${ann_vcf}.gz" into ann_vcf_zip_ch
    path "${ann_vcf}.gz.tbi" into ann_vcf_ch_index

    script:
    """
    bgzip $ann_vcf
    tabix -p vcf ${ann_vcf}.gz
    """
}

//Combine SNP and Indel VCFs
process combine_variants {
    publishDir params.outdir, mode:'copy'

    input:
    path ann_vcf from ann_vcf_zip_ch.collect()
    path ann_vcf_index from ann_vcf_ch_index.collect()

    output:
    path "combined_filtered_pass.vcf.snpeff.vcf.gz" into comb_ann_vcf_ch
    
    script:
    """
    bcftools concat -a -o combined_filtered_pass.vcf.snpeff.vcf.gz -O z ${ann_vcf.join(' ')}
    """
}
