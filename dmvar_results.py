#!/usr/bin/env python3

# dmvar_results.py

# AUTHOR: Adam James Reid
# Copyright (C) 2022 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

# NOTES 
# Better written version of dmvar_results.py
 # Make genesummary data optional
 # Allow gene positional information to be added from file
 # Potential for identifying overlaps with deficiency lines

# N.b. the funciton filter currently works in only a limited way - if not specified there is no filter. HIGH and MODERATE can be specified to identify, by priority HIGH effect annotations for variantsand if absent, MODERATE effect annotations e.g. using '-f HIGH,MODERATE'. Or just '-f HIGH' for only HIGH effect variants. If specified, variantswith neither of these annotations will be ignored.

# imports
import sys
import subprocess
#import ExcelWriter, DataFrame from pandas
from pandas import ExcelWriter, DataFrame
import openpyxl
import argparse
import re
import gzip

###########
# FUNCTIONS

def parse_samplesheet(ss_file):
    '''Return dictionary of lists of parent->progeny and dictionary of tuples of progent ->het/hom, chromosome of interest'''
    progeny = dict()
    details = dict()

    with open(ss_file) as ss:
        for x in ss.readlines():
            x = x.rstrip()
            v = x.split(',')

            # Ignore header line
            if v[0] == 'sample':
                continue

            # If v[3] is 'None', this is a parent, we can get all info from progeny lines 
            if v[3] != 'None':
                if v[3] not in progeny:
                    progeny[v[3]] = list()
                # Add sample name v[0] under parent v[3]
                progeny[v[3]].append(v[0])

                if v[0] not in details:
                    details[v[0]] = tuple()
                # Add het/hom and chromosome of interest under sample v[0]
                details[v[0]] = (v[1], v[2])

    return progeny, details

def read_file(fn):
    '''Generator for reading files either zipped or not'''
    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            for line in f:
                yield line
    else:
        with open(fn) as f:
            for line in f:
                yield line

def print_vcf_header(vcf, fh):
    ''' Write VCF header out to file specified by fh'''
    for x in read_file(vcf):
        if x.startswith('#'):
            x = x.rstrip()
            fh.write('{}\n'.format(x))

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

def is_hom_fuzzy(gt, af, min_hom_allele_freq):
    ''' Is this genotype homozygous? It can have been called as het, but as long as the frequency of the alternative allele is >= min_hom_allele_freq we will call it homozygous anyway '''
    # If min_hom_allele_freq == False, do a normal is_hom call
    if min_hom_allele_freq == False:
        print("Run is_hom")
        return is_hom(gt)
    # If it is a het, check allele frequencies, return true if alternative allele frequency is high enough
    elif is_het(gt):
        allele_freqs = af.split(',')
        # Determine total depth over this variant to check that it is not zero
        depth_sum = sum([float(x) for x in allele_freqs])
        
        if depth_sum > 0:
            # frequency of alternate allele, assuming a single alternate
            alt_freq = float(allele_freqs[1]) / (sum([float(x) for x in allele_freqs]))
        else:
            return False

        print ("{} {} {} {}".format(gt, af, alt_freq, min_hom_allele_freq))
        if alt_freq >= min_hom_allele_freq:
            print ("True")
            return True
        else:
            return False
    # Otherwise return False
    else:
        return False

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

def single_progeny_het_fuzzy(sample_calls, parent, max_shared):
    ''' Determine whether a subset of progeny are heterozygous for a particular variant - alternative to single_progeny_het'''
    het_progeny = list()
    for s in sample_calls:
        if s != parent:
            if(is_het(sample_calls[s])):
                het_progeny.append(s)
    if len(het_progeny) <= (int(max_shared) + 1):
        return het_progeny
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

def single_progeny_hom(sample_calls, parent):
    ''' Determine whether a single progeny is homozygous alt or ref, while parent and other progeny are homozygous but opposite allele '''
    # Determine if parent is homozygous and alt or ref
    p_alleles = re.split('[|/]', sample_calls[parent])
    p_gt = '0'
    if p_alleles[0] == '0' and p_alleles[1] == '0':
        p_gt = '0'
    elif p_alleles[0] == '1' and p_alleles[1] == '1':
        p_gt = '1'
    # Return None if parent is heterozygous
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

def single_progeny_hom_fuzzy(sample_calls, sample_freqs, parent, max_shared, min_hom_allele_freq):
    ''' Alternative funciton to single_progeny_hom which allows heterozygous calls with alternative allele freq >= min_allele_freq 
    to be called as homozygous
    '''

    print("Running single_progeny_hom_fuzzy {} {} {}".format(parent, max_shared, min_hom_allele_freq))

    # Determine if parent is homozygous and alt or ref
    p_alleles = re.split('[|/]', sample_calls[parent])
    p_gt = '0'
    if p_alleles[0] == '0' and p_alleles[1] == '0':
        p_gt = '0'
    elif p_alleles[0] == '1' and p_alleles[1] == '1':
        p_gt = '1'
    # Return None if parent is heterozygous
    else:
        print("Parent is heterozygous")
        return None

    res_dict = {"hom_p" : 0, "hom_np" : 0, "het" : 0}
    hom_np_samples = list()

    # Go through calls 
    for s in sample_calls:
        if s != parent:
            s_alleles = re.split('[|/]', sample_calls[s])
            # If this is a heterozygous call
            if s_alleles[0] != s_alleles[1] and '.' not in s_alleles: # Het call
                print("Looks like a het {} {}".format(s_alleles[0], s_alleles[1]))
                # if this is a fuzzy homozygote
                if is_hom_fuzzy(sample_calls[s], sample_freqs[s], min_hom_allele_freq):
                    print("is_hom_fuzzy")
                    # if the fuzzy homozygote is not the parental genotype
                    if s_alleles[1] != p_gt:
                        res_dict['hom_np'] = res_dict['hom_np'] + 1
                        hom_np_samples.append(s)
                    # It is a homozygous parental genotype
                    else:
                        res_dict['hom_p'] = res_dict['hom_p'] + 1
                # It is not a fuzzy hom, so it is a het
                else:
                    print("not a fuzzy hom, but a het")
                    res_dict['het'] = res_dict['het'] + 1
            # Homozygous, parental genotype
            elif s_alleles[0] == s_alleles[1] and '.' not in s_alleles and s_alleles[0] == p_gt: # Hom parental call
                res_dict['hom_p'] = res_dict['hom_p'] + 1
            # Homozygous alternative
            elif s_alleles[0] == s_alleles[1] and '.' not in s_alleles and s_alleles[0] != p_gt: # Hom alternate call
                res_dict['hom_np'] = res_dict['hom_np'] + 1
                hom_np_samples.append(s)
    
    # If the number of homozygous nonparentals is gt 0 (i.e. we have at least one interesting variant and the number of hets and hom_nps is not more than the max_shared+1, return them 
    print("res_dict {}".format(res_dict))
    if res_dict['hom_np'] > 0 and ((res_dict['hom_np'] + res_dict['het']) <= (max_shared + 1)):
        print("Return 1 : {}".format(hom_np_samples))
        return hom_np_samples
    else:
        print("Return 0 : {}".format(hom_np_samples))
        return None


def parse_function(info):
    ''' Given an info line, presumably tagged by snpeff, identify ANN records and return gene name, mutation type etc'''
    ''' 'effects' can contain (HIGH, MODERATE, LOW, MODIFIER)'''
    return_data = dict()
    fields = info.split(';')

    for f in fields:
        if f.startswith('ANN='):
            f = f.replace('ANN=', '')
            anns = f.split(',')
            anns_join = "\n".join(anns)

            for a in anns:
                abits = a.split('|')
                if abits[2] not in return_data:
                    return_data[abits[2]] = list()
                index_names = [0,1,2,3,4,9,10]
                return_data[abits[2]].append([abits[val] for val in index_names])

    return return_data

def read_gene_summary(fn):
    ''' Read gene summary data into a dict '''
    info = dict()
    #with open(fn, 'rt') as f:
    for x in read_file(fn):
        if x.startswith('#'):
            continue
        v = x.split('\t')
        v[1].rstrip() # Sometimes there seem to be trailing newlines
        info[v[0]] = v[1]
    return info

def read_genemap(fn):
    '''Read gene map file from Flybase into a dict'''
    info = dict()
    for x in read_file(fn):
        if x.startswith('#') or x == '':
            continue
        x = x.rstrip()
        v = x.split('\t')
        # sixth column in the sequence locus
        # third column is the primary_FBid
        if len(v) > 5:
            info[v[2]] = v[5]
    return info

def function_effect(info):
    '''
    Determine the most serious kind of effect the variant has
    Return the seriousness of the effect and the details
    '''
    # Use parse_function list to get a dictionary of effect type to lists of effect details
    func_vals = parse_function(info)

    # Return is a list of effect summary data
    out_list = list()

    if func_vals:
        effect_type = ''
        if 'HIGH' in func_vals:
            effect_type = 'HIGH'
        elif 'MODERATE' in func_vals:
            effect_type = 'MODERATE'
        elif 'LOW' in func_vals:
            effect_type = 'LOW'
        elif 'MODIFIER' in func_vals:
            effect_type = 'MODIFIER'
        else:
            sys.stderr.write('No effect type. SnpEff annotation missing? {}'.format(info))

        out_list.append('\t'.join(func_vals[effect_type][0]))
            
        # Determine gene name if any
        gene = func_vals[effect_type][0][4]
        if gene:
            out_list.append(gene)
        else:
            out_list.append("")

        return (effect_type, out_list)
    else:
        return (None, None)
    

def select_results(f, prog, chrom):
    '''
    Select relevant results for a particular mutant from the annotated variants results file
    f = file of annotated results, either hom or het variants
    prog = progeny/sample name
    chrom = chromosome/chromosome arm of interest
    '''
    results = list()
    with open(f) as fh:
        for x in fh.readlines():
            x = x.rstrip()
            v = x.split('\t')
            # filter for the progeny name and chromosome of interest.
            if v[0] == prog and v[1].startswith(chrom):
                results.append(v)
    return results


###########
# argparse
parser = argparse.ArgumentParser(description='Idenitify variants of interest from Drosophila mutatation experiments')
parser.add_argument('-g', '--genesummary', help="File downloaded from Flybase with detailed info about each gene [optional]")
parser.add_argument('-m', '--genemapfile', help="Gene map file from Flybase with locations of genes [optional]")
parser.add_argument('-s', '--samplesheet', help="Samplesheet with header 'sample,chromosome,type,control,batch,vcf'")
parser.add_argument('-v', '--vcf_file', help="VCF file of variants from dmvar.nf pipeline - including all samples in samplesheet")
parser.add_argument('-f', '--funcfilter', help="Include only SNPs with annotations labelled as one or more of HIGH, MODERATE, LOW, MODIFIER - comma separated list")
parser.add_argument('--maxshared', help="How many siblings can share an alternative allele and still be called [0]")
parser.add_argument('-a', '--minhomfreq', help="Turn on fuzzy homozygous calling and set minimum alternative allele frequency to call a homozygote [None]")
args = parser.parse_args()

###################
# manage arguments
if args.samplesheet:
    ss_file = args.samplesheet
else:
    print("Please provide a samplesheet (-s)")
    exit()
if args.genesummary:
    #genesummary_file = args.genesummary
    genesummary = read_gene_summary(args.genesummary)
if args.genemapfile:
    genemap = read_genemap(args.genemapfile)
if args.vcf_file:
    vcf_file = args.vcf_file
else:
    print("Please provide a VCF file of variants from dmvar.nf pipeline (-v)")
    exit()
if args.funcfilter:
    funcfilter = args.funcfilter.split(',')
if args.maxshared:
    max_shared = int(args.maxshared)
else:
    max_shared = 0
if args.minhomfreq:
    min_hom_allele_freq = float(args.minhomfreq)
else:
    min_hom_allele_freq = False



#################
# Main

# Dict of progeny for each parent, Dict of het/hom and chromosome of interest for each progeny
progeny, details = parse_samplesheet(ss_file)

# for each parent
for parent in progeny:

    # define cohort subset vcf
    cohort_vcf_file = vcf_file+'.'+parent+'.vcf'

    # Define cohort homozygous filtered vcf
    cohort_vcf_file_hom = vcf_file+'.'+parent+'_hom.vcf'
    homout_vcf = open(cohort_vcf_file_hom, 'w')
    print_vcf_header(cohort_vcf_file_hom, homout_vcf)

    # Define cohort heterozygous filtered vcf
    cohort_vcf_file_het = vcf_file+'.'+parent+'_het.vcf'
    hetout_vcf = open(cohort_vcf_file_het, 'w')
    print_vcf_header(cohort_vcf_file_het, hetout_vcf)

    # Define homozygous results file
    cohort_hom_variants_file = vcf_file+'.'+parent+'.hom.txt'
    homout = open(cohort_hom_variants_file, 'w')
    homout.write("Sample\tChromosome\tPosition\tAlt allele\tMutation type\tMutation effect\tGene\tFB gene id\tNT change\tAA change\tShared_with\tGene location\tGene description\n")

    # Define heterozygous results file
    cohort_het_variants_file = vcf_file+'.'+parent+'.het.txt'
    hetout = open(cohort_het_variants_file, 'w')
    hetout.write("Sample\tChromosome\tPosition\tAlt allele\tMutation type\tMutation effect\tGene\tFB gene id\tNT change\tAA change\tShared with\tGene location\tGene description\n")
    
    # bcftools to subset vcf to cohort
    cmd = "which bcftools; bcftools -v;bcftools view -s {},{} {} > {}".format(','.join(progeny[parent]), parent, vcf_file, cohort_vcf_file)
    print(cmd)
    returned_value = subprocess.call(cmd, shell=True)

    # identify hom/het snps of interest
    # Get list of samples from the cohort VCF file
    samples = get_vcf_sample_names(cohort_vcf_file)

    # read cohort vcf
    for x in read_vcf(cohort_vcf_file):
        chrom = x[0]
        pos = x[1]

        # Get genotypes for each sample
        genotypes = [s.split(':')[0] for s in x[9:]]

        # Get allele_freqs
        allele_freqs = [s.split(':')[1] for s in x[9:]]

        # Make a dict of sample->genotype
        sample_calls = dict(zip(samples, genotypes))
        # Make a dict of sample->allele freq
        sample_freqs = dict(zip(samples, allele_freqs))

        # If parent is homozygous (otherwise we are not interested)
        if is_hom(sample_calls[parent]):
            # is there a single (or max_shared) progeny with a heterozygous non-parental call?
            p = single_progeny_het_fuzzy(sample_calls, parent, max_shared)
            if p:
                out_line = ''
                genesummary_text = ''
                genelocus_text = ''
                # Alternative here to parse function effect first, see what we get back and whether that is in the filter
                #if args.funcfilter:
                # parse function - pass annotatation field from VCF and comma-sep list of effects accepted e.g. 'HIGH,MODERATE'
                (effect_type, effect) = function_effect(x[7])
                # Continue with the variant if the most serious effect type found is in the supplied list, or if there is no funcfilter
                if not args.funcfilter or effect_type in funcfilter:
                    if effect[1] and args.genesummary and effect[1] in genesummary:
                        genesummary_text = genesummary[effect[1]]
                    if effect[1] and args.genemapfile and effect[1] in genemap:
                        genelocus_text = genemap[effect[1]]

                    # Determine how many progeny share this variant
                    if len(p) > 1:
                        shared_with = len(p) - 1
                    else:
                        shared_with = 0

                    # Loop through return progeny and write out variants
                    for i in p:
                        out_line = '\t'.join([i, chrom, pos, effect[0], str(shared_with), genelocus_text, genesummary_text])
                        hetout.write('{}\n'.format(out_line))

                # write out VCF line
                vcf_line = "\t".join(x)
                hetout_vcf.write('{}\n'.format(vcf_line))

            # is there a single (or max_shared) progeny with a homozygous (or fuzzy homozygous) non-parental call?
            #print("Variant site: {} {}".format(chrom, pos))
            p = single_progeny_hom_fuzzy(sample_calls, sample_freqs, parent, max_shared, min_hom_allele_freq)
            if p:
                print("{} = p from single_progeny_hom_fuzzy".format(p))
                out_line = ''
                genesummary_text = ''
                genelocus_text = ''
                # parse function (n.b. effect[1] is the gene name)
                (effect_type, effect) = function_effect(x[7])
                #if effect_type in funcfilter or not args.funcfilter:
                if not args.funcfilter or effect_type in funcfilter:
                    if effect[1] and args.genesummary and effect[1] in genesummary:
                        genesummary_text = genesummary[effect[1]]
                    if effect[1] and args.genemapfile and effect[1] in genemap:
                        genelocus_text = genemap[effect[1]]

                    # Determine how many progeny share this variant
                    if len(p) > 1:
                        shared_with = len(p) - 1
                    else:
                        shared_with = 0

                    # Loop through return progeny and write out variants
                    for i in p:
                        out_line = '\t'.join([i, chrom, pos, effect[0], str(shared_with), genelocus_text, genesummary_text])
                        homout.write('{}\n'.format(out_line))

                # write out VCF line
                vcf_line = "\t".join(x)
                homout_vcf.write('{}\n'.format(vcf_line))
    
    # Get gene position information from file if present

    # Close intermediate cohort-specific results files
    homout_vcf.close()
    hetout_vcf.close()
    hetout.close()
    homout.close()
            
    # Initialise spreadsheet for this cohort
    excel_file = parent + '.xlsx'
    writer = ExcelWriter(excel_file)

    # for each progeny - get relevant variants from file
    for prog in progeny[parent]:
        # get lines of file relating to relevant progeny, het/hom and chromosome
        if details[prog][1] == 'Heterozygous':
            prog_res_list = select_results(cohort_het_variants_file, prog, details[prog][0])
        elif details[prog][1] == 'Homozygous':
            prog_res_list = select_results(cohort_hom_variants_file, prog, details[prog][0])
        else:
            print("Not a valid analysis type, should be Heterozygous or Homogygous: {} {}".format(prog, details[prog][1]))
            exit()
        
        # Write out Excel sheet
        #if args.funcfilter:
        df = DataFrame(prog_res_list, columns=['Sample','Chromosome','Position', 'Alt allele', 'Mutation type', 'Mutation effect', 'Gene', 'FB gene id', 'NT change', 'AA change', 'Shared with', 'Gene location', 'Gene description'])
        #else:
        #    df = DataFrame(prog_res_list, columns=['Sample','Chromosome','Position'])
        sheet_name = prog+'_'+details[prog][0]+'_'+details[prog][1]
        df.to_excel(writer, index=False, sheet_name=sheet_name)

        # Output each mutant to a tsv file as well as to Excel
        df.to_csv(sheet_name+'.tsv', sep='\t', index=False)

    writer.save()
