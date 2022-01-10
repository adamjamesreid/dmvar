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
import argparse

parser = argparse.ArgumentParser(description='Idenitify variants of interest from Drosophila mutatation experiments')
parser.add_argument('-g', '--genesummary', help="File downloaded from Flybase with detailed info about each gene")
parser.add_argument('-s', '--samplesheet', help="Samplesheet with header 'sample,chromosome,type,control,batch,vcf'")
parser.add_argument('-r', '--remove_parental_script', help="Path to remove_parental2.py script")
parser.add_argument('-v', '--vcf_file', help="VCF file of variants from dmvar.nf pipeline - including all samples in samplesheet")
args = parser.parse_args()

if args.samplesheet:
    ss_file = args.samplesheet
else:
    print("Please provide a samplesheet (-s)")
    exit()
if args.genesummary:
    genesummary_file = args.genesummary
if args.remove_parental_script:
    remove_parental_script = args.remove_parental_script
else:
    print("Please provide a path to the remove_parental2.py script (-r)")
    exit()
if args.vcf_file:
    vcf_file = args.vcf_file
else:
    print("Please provide a VCF file of variants from dmvar.nf pipeline (-v)")
    exit()


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

for parent in progeny:
    cohort_vcf_file = vcf_file+'.'+parent+'.vcf'
    cohort_hom_variants_file = vcf_file+'.'+parent+'.homs.txt'
    cohort_het_variants_file = vcf_file+'.'+parent+'.hets.txt'

    cmd = "bcftools view -s {},{} {} > {}".format(','.join(progeny[parent]), parent, vcf_file, cohort_vcf_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

    cmd = "python ~/code_development/dmvar/remove_parental2.py -p {} -o {}_cohortfilt_snps_homs.vcf -a 2 -g {} -f HIGH,MODERATE {} > {}".format(parent, parent, genesummary_file, cohort_vcf_file, cohort_hom_variants_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

    cmd = "python ~/code_development/dmvar/remove_parental2.py -p {} -o {}_cohortfilt_snps_homs.vcf -a 1 -g {} -f HIGH,MODERATE {} > {}".format(parent, parent, genesummary_file, cohort_vcf_file, cohort_het_variants_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

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
            print("Not a valid analysis type, should be Heterozygous or Homogygous: {} {}".format(prog, details[prog][1]))
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
