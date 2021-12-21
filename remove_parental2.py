# Take a processed vcf of just a single parental lineage and identify SNPs which are present in the parent and filter these out

# Can call either:
    # 1. Variants which are homozygous in parent and heterozygous in a single progeny
    # 2. Homozygous in the parent and homozygous but non-parental in only a single progeny. All other progeny have homozygous parental genotype
    # n.b. uncalled alleles are ignored

# Functional annotation:
    # Using Snpeff annotations, output variants with HIGH impact and/or NMD and/or LOF 

#####
# To do in version 2
########

# 1. Remove NMD/LOF and just use ANN
# 2. Get gene name from ANN (assume just one gene?) and use this for getting FlyBase info
# 3. Output nice Excel file

import sys
import gzip
import re
import argparse

'''
Generator for reading files either zipped or not
'''
def read_file(fn):
    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            for line in f:
                yield line
    else:
        with open(fn) as f:
            for line in f:
                yield line

def get_vcf_sample_names (vcf):
    ''' Return a list of sample names for a vcf '''
    for x in read_file(vcf):
        if x.startswith('#CHROM'):
            x = x.rstrip()
            v = x.split('\t')
            samples = v[9:]

            return samples

def read_vcf (vcf):
    ''' Generator to read each non-header line of a vcf in as a list of elements'''
    for x in read_file(vcf):
        if x.startswith('#'):
            continue
        else:
            x = x.rstrip()
            v = x.split('\t')
            yield(v)

def is_hom (gt):
    ''' Is this genotype homozygous? If any allele is uncalled "." then it is not homozygous'''
    alleles = re.split('[|/]', gt)
    # Check whether all genotype values are the same
    allele_set = set()
    for a in alleles:
        allele_set.add(a)
    if len(allele_set) == 1 and '.' not in allele_set:
        return True
    else:
        return False

# Only just started this function
def genotypes_chromosomes (vcf):
    ''' Determine the numbers of different genotypes per sample, per chromosome '''
    samples = get_vcf_sample_names(vcf)

    for x in read_vcf(vcf):
        chrom = x[0]
        sample_calls = dict(zip(samples, x[9:]))
        print(sample_calls)
        

def single_progeny_hom(sample_calls, parent):
    ''' Determine whether a single progeny is homozygous alt or ref, while parent and other progeny are homozygous but opposite allele '''
    # Determine if parent is homozygous and alt or ref
    p_alleles = re.split('[|/]', sample_calls[parent])
    p_gt = '0'
    if p_alleles[0] == '0' and p_alleles[1] == '0':
        p_gt = '0'
    elif p_alleles[0] == '1' and p_alleles[1] == '1':
        p_gt = '1'
    else:
        return None

    res_dict = {"hom_p" : 0, "hom_np" : 0, "het" : 0}
    hom_np_samples = list()

    # Determine progeny alleles
    for s in sample_calls:
        if s != parent:
            s_alleles = re.split('[|/]', sample_calls[s])
            if s_alleles[0] != s_alleles[1] and '.' not in s_alleles: # Het call
                res_dict['het'] = res_dict['het'] + 1
            elif s_alleles[0] == s_alleles[1] and '.' not in s_alleles and s_alleles[0] == p_gt: # Hom parental call
                res_dict['hom_p'] = res_dict['hom_p'] + 1
            elif s_alleles[0] == s_alleles[1] and '.' not in s_alleles and s_alleles[0] != p_gt: # Hom alternate call
                res_dict['hom_np'] = res_dict['hom_np'] + 1
                hom_np_samples.append(s)

    if res_dict['het'] == 0 and res_dict['hom_np'] == 1:
        return hom_np_samples[0]
    else:
        return None

def single_progeny_het(sample_calls, parent):
    ''' Determine whether a single progeny is heterozygous, while others are homogygous '''
    hets = 0
    progeny_to_return = ''
    for s in sample_calls:
        if s != parent:
            if(is_het(sample_calls[s])):
                hets = hets + 1
                progeny_to_return = s # this will be returned if only one progreny is a het
    if hets == 1:
        return progeny_to_return
    else:
        return None

def is_het(gt):
    alleles = re.split('[|/]', gt)
    # Check whether all genotype values are the same
    allele_set = set()
    for a in alleles:
        allele_set.add(a)
    ###!!!! N.b. for multiallelic sites this might fail!!!!###
    if len(allele_set) > 1 and '.' not in allele_set:
        return True
    else:
        return False

def print_vcf_header(vcf, fh):
    ''' Write VCF header out to file specified by fh'''
    for x in read_file(vcf):
        if x.startswith('#'):
            x = x.rstrip()
            fh.write('{}\n'.format(x))

def parse_function(info, effects):
    ''' Given an info line, presumably tagged by snpeff, identify ANN records and return gene name, mutation type etc''' 
    ''' 'effects' can contain (HIGH, MODERATE, LOW, MODIFIER)'''
    return_data = dict()
    fields = info.split(';')
    #print(fields, "\n")
    for f in fields:
        if f.startswith('ANN='):
            f = f.replace('ANN=', '')
            anns = f.split(',')
            anns_join = "\n".join(anns)
            #print(anns_join, "\n")
            for a in anns:
                abits = a.split('|')
                if abits[2] not in return_data:
                    return_data[abits[2]] = list()
                index_names = [0,1,2,3,4,9,10]
                return_data[abits[2]].append([abits[val] for val in index_names])
    #print(return_data,"\n")
    return return_data

def read_gene_summary(fn):
    ''' Read gene summary data into a dict '''
    info = dict()
    with open(fn, 'rt') as f:
        for x in f.readlines():
            if x.startswith('#'):
                continue
            v = x.split('\t')
            v[1].rstrip() # Sometimes there seem to be trailing newlines 
            info[v[0]] = v[1]
    return info


# Parse command line arguments #

parser = argparse.ArgumentParser(description='Filter a Drosophila cohort VCF')
parser.add_argument('-p', '--parent', help="Sample name for parent")
parser.add_argument('vcf', help="VCF file containing genotyped/filtered SNPs for cohort, can be bgzipped")
parser.add_argument('-o', '--outfile', help="Name of VCF output file [cohort_filt.vcf]", default="cohort_filt.vcf")
parser.add_argument('-f', '--funcfilter', help="Include only SNPs with annotations labelled as one or more of HIGH, MODERATE, LOW, MODIFIER - comma separated list")
parser.add_argument('-a', '--allelefilter', type=int, help="Include only alleles which are homozygous in parent and heterozygous in a single progeny (1) or homozygous in the parent and homozygous but different from the parent in only a single progeny (2) default = [1]")
parser.add_argument('-g', '--genesummary', help="File downloaded from Flybase with detailed info about each gene")
args = parser.parse_args()

if args.parent:
    parent = args.parent
if args.vcf:
    vcf = args.vcf
if args.outfile:
    outfile = args.outfile
if args.genesummary:
    # Read gene summary data
    genesummary = read_gene_summary(args.genesummary)
else:
    genesummary = dict()
if args.funcfilter:
    funcfilter = args.funcfilter.split(',')
if args.allelefilter:
    allele_filter = args.allelefilter
else:
    allele_filter = 1
            
# possible outputs
## VCF filtered for SNPs of interest e.g. those in hom in parent and het in a single progeny
## VCF filtered for SNPs of interest e.g. those in hom in parent and het in a single progeny and LOF/NMD
## Table of:
    # Sample name
    # Chromosome arm
    # Mutation details
    # Gene name
    # Gene description

# Setup VCF output
outh = open(outfile, 'w')
print_vcf_header(vcf, outh)

# Print header for STDOUT
print("Sample\tChromosome\tPosition\tAlt allele\tMutation type\tMutation effect\tGene\tFB gene id\tNT change\tAA change\tGene description")

samples = get_vcf_sample_names(vcf)

for x in read_vcf(vcf):
    chrom = x[0]
    pos = x[1]
    # Get genotypes for each sample
    genotypes = [s.split(':')[0] for s in x[9:]]

    sample_calls = dict(zip(samples, genotypes))

    # Parent is homozygous
    if is_hom(sample_calls[parent]):
        r = None
        # one sibling only is heterozygous
        if allele_filter == 1:
            r = single_progeny_het(sample_calls, parent)
        # one sibling only is homozygous non-parental, while parent is homozygous, no hets
        elif allele_filter == 2:
            r = single_progeny_hom(sample_calls, parent)

        if r:
            out_list = [r, chrom, pos]

            # Functional filtering (e.g. look for HIGH/MODERATE effects)
            # n.b. we return just the first annotation, preferntially with HIGH effect
            if args.funcfilter:
                func_vals = parse_function(x[7], funcfilter)
                if func_vals:
                    if 'HIGH' in func_vals:
                        out_list.append('\t'.join(func_vals['HIGH'][0]))

                        # Determine gene name if any
                        gene = func_vals['HIGH'][0][4]
                        if gene and gene in genesummary:
                            out_list.append(genesummary[gene].strip()) # N.b. there are sometimes trailing newlines
                        else:
                            out_list.append("")
                        
                        print('\t'.join(out_list))
                    elif 'MODERATE' in func_vals:
                        out_list.append('\t'.join(func_vals['MODERATE'][0]))

                        # Determine gene name if any and add Flybase gene summary
                        gene = func_vals['MODERATE'][0][4]
                        if gene and gene in genesummary:
                            out_list.append(genesummary[gene].rstrip()) # N.b. there are sometimes trailing newlines
                        else:
                            out_list.append("")

                        print('\t'.join(out_list))

                    vcf_line = "\t".join(x)
                    outh.write('{}\n'.format(vcf_line))
            else:
                vcf_line = "\t".join(x)
                outh.write('{}\n'.format(vcf_line))

                print('\t'.join(out_list))

        
outh.close()
