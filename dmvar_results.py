# Take samplesheet, identify each parent (none in parent column), subset other samples per parent
# subset VCFs based on cohorts
# Determine whether to do hets, homs or both
# Run remove_parental2.py using hets/homs
# For each samples, filter variants based on sample, het/hom, chromosome of interest
# Write out as a spreadsheet

import sys
import subprocess
import pandas as pd
import openpyxl

ss_file = sys.argv[1]
vcf_file = sys.argv[2]

# list of parents
# dict of parent -> progeny
progeny = dict()
# dict of progeny -> het/hom, chromosome of interest
details = dict()

with open(ss_file) as ss:
    for x in ss.readlines():
        x = x.rstrip()
        v = x.split(',')
        
        if v[0] == 'sample':
            continue
        
        if v[3] != 'None':
            if v[3] not in progeny:
                progeny[v[3]] = list()
            progeny[v[3]].append(v[0])
            
            if v[0] not in details:
                details[v[0]] = tuple()
            details[v[0]] = (v[1], v[2])

#print(progeny)

#print(details)

for parent in progeny:
    cohort_vcf_file = vcf_file+'.'+parent+'.vcf'
    cohort_hom_variants_file = vcf_file+'.'+parent+'.homs.txt'
    cohort_het_variants_file = vcf_file+'.'+parent+'.hets.txt'

    cmd = "bcftools view -s {},{} {} > {}".format(','.join(progeny[parent]), parent, vcf_file, cohort_vcf_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

    cmd = "python ~/code_development/dmvar/remove_parental2.py -p {} -o {}_cohortfilt_snps_homs.vcf -a 2 -g ~/code_development/dmvar/automated_gene_summaries.tsv -f HIGH,MODERATE {} > {}".format(parent, parent, cohort_vcf_file, cohort_hom_variants_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

    cmd = "python ~/code_development/dmvar/remove_parental2.py -p {} -o {}_cohortfilt_snps_homs.vcf -a 1 -g ~/code_development/dmvar/automated_gene_summaries.tsv -f HIGH,MODERATE {} > {}".format(parent, parent, cohort_vcf_file, cohort_het_variants_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)
    #python ~/code_development/dmvar/remove_parental2.py -p Actrl -o Actrl_cohortfilt_snps_homs.vcf -a 2 -g ~/code_development/dmvar/automated_gene_summaries.tsv -f HIGH,MODERATE snps_filtered_pass.vcf.snpeff.Actrl.vcf  > Actrl_cohortfilt_snps_homs.txt

    # Declare a spreadsheet output for this cohort
    excel_file = parent + '.xlsx'
    #writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
    writer = pd.ExcelWriter(excel_file)
    #writer.book = openpyxl.load_workbook(excel_file)

    for prog in progeny[parent]:
        prog_res_list = list()
        # pick het/hom file
        file_to_open = ''
        if details[prog][1] == 'Heterozygous':
            file_to_open = cohort_het_variants_file
        elif details[prog][1] == 'Homozygous':
            file_to_open = cohort_hom_variants_file
        else:
            print("Not a valid analysis type, should be Heterozygous or Homogygous: {}".format(details[prog][1]))
            exit()

        with open(file_to_open) as hethom_file:
            for x in hethom_file.readlines():
                x = x.rstrip()
                v = x.split('\t')
                # filter for the progeny name and chromosome of interest.
                if v[0] == prog and v[1].startswith(details[prog][0]):
                    prog_res_list.append(v)
        
        # Output to an excel file
        df = pd.DataFrame(prog_res_list, columns=['Sample','Chromosome','Position', 'Alt allele', 'Mutation type', 'Mutation effect', 'Gene', 'FB gene id', 'NT change', 'AA change', 'Gene description'])
        #print (prog_res_list)
        #print(df)
        df.to_excel(writer, index=False, sheet_name=prog+'_'+details[prog][0]+'_'+details[prog][1])
    writer.save()