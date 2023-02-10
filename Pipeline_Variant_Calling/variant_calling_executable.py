#### Imports ####

import pybedtools as pb
import subprocess as sp
import os
from pathlib import Path
import argparse
import time
import pandas as pd
import numpy as np
from itertools import chain
from collections import defaultdict

### FUNTIONS ####

def index_feature_file(gatk, feat_file):
    cmd = '{} IndexFeatureFile -I {}' .format(gatk, feat_file)
    sp.call(cmd, shell = True)

def AddOrReplaceReadGroups(gatk, bam, outfld):
    outfile = Path(outfld).joinpath(Path(bam).stem+'_withRG.bam')
    name = Path(bam).stem
    cmd = ("{} "
           "AddOrReplaceReadGroups "
           "-INPUT {} "
           "-OUTPUT {} "
           "-RGID group_{} "
           "-RGLB lib_{} "
           "-RGPL illumina "
           "-RGPU unit1 "
           "-RGSM {}_sample") .format(gatk, bam, outfile, name, name, name)
    sp.call(cmd, shell=True)

def mark_duplicates(gatk, bam, outfld):
    outfile = Path(outfld).joinpath(Path(bam).stem+'_markedDuplicates.bam')
    mtrcsfile = Path(outfld).joinpath(Path(bam).stem+'_metrics.txt')
    args = [gatk, bam, outfile, mtrcsfile]
    cmd = '{} MarkDuplicates -I {} -O {} -M {}' .format(*args)
    sp.call(cmd, shell=True)

    cmd = 'samtools sort {} -o {}' .format(outfile, outfile)
    sp.call(cmd, shell=True)

def base_recalibration(gatk, bam, outfld, known, ref_fasta):
    outfile = Path(outfld).joinpath(Path(bam).stem+'_baserecal_table.table')
    args = [gatk, bam, ref_fasta, known, outfile]
    cmd = ('{} BaseRecalibrator '
           '-I {} -R {} '
           '--known-sites {} '
           '-O {}') .format(*args)
    sp.call(cmd, shell=True)

def applyBQSR(gatk, bam, outfld, ref_fa):

    outfile = Path(outfld).joinpath(Path(bam).stem+'_BQSR.bam')
    recal_table = Path(outfld).joinpath(Path(bam).stem+'_baserecal_table.table')
    args = [gatk, bam, ref_fa, recal_table, outfile]

    cmd = ('{} ApplyBQSR '
           '-I {} -R {} '
           '--bqsr-recal-file {} '
           '-O {}') .format(*args)

    sp.call(cmd, shell=True)

def mergeBams(gatk, *bams, outfile):
    nbams = len(bams)
    inputs = '-I {} '*nbams
    cmd = '{} ' .format(gatk) + \
        'MergeSamFiles ' + \
        inputs .format(*bams) + \
        '-O {}.bam' .format(outfile)
    sp.call(cmd, shell=True)

def samtools_index(samtools_path, infile):
    cmd = [samtools_path, 'index', infile]
    sp.run(cmd)

def call_variants(gatk, ref_fasta, bam, outfile):
    args = [gatk, ref_fasta, bam, outfile]
    cmd = ('{} HaplotypeCaller '
           '-R {} -I {} -O {} -ploidy 1') .format(*args)

    sp.call(cmd, shell=True)

def call_VEP(vep, vcf, gff, fasta):

    outfile = Path(vcf).parent.joinpath(Path(vcf).stem+'_VEPannotated.txt')
    args = [vep, vcf, outfile, gff, fasta]

    cmd = ("{} "
           "-i {} "
           "-o {} "
           "--gff {} "
           "--fasta {} "
           "--force_overwrite "
           "--vcf") .format(*args)

    sp.call(cmd, shell=True)

def getAnnot(gffentry):

    info = gffentry.fields[8].split(";")
    dinfo = {x.split('=')[0]:x.split('=')[1] for x in info}
    gid = dinfo['ID']
    anot = dinfo['description']
    return([gid, anot])

def getRatioDepth(GF):
    if len(GF) <2:
        rf = np.nan
        alt = np.nan
        ratio = np.nan
        dp = 0
    else:
        rf = int(GF[1].split(",")[0])
        alt = int(GF[1].split(",")[1])
        dp = rf+alt

        if dp == 0:
            ratio = np.nan
        else:
            ratio = round(rf / dp, 2)

    return(rf, alt, ratio, dp)

def parse_variant(variant, annot_dict):

    # Parse vcf info
    ref = variant.fields[3]
    alt = variant.fields[4]
    pos = variant.start
    chrom = variant.chrom

    samples = variant.fields[9:]
    sample_vars = [x.split(':') for x in samples]

    col_vals = [getRatioDepth(x) for x in sample_vars]
    unnested_col_vals = [item for sublist in col_vals for item in sublist]
    #ref_count1, alt_count1, r1, d1 = getRatioDepth(v10G)

    parsed_vcf = [chrom, pos, ref, alt]+unnested_col_vals

    # Parse vep info
    info = {}
    for x in variant.fields[7].split(";"):
        feat = x.split("=")
        if len(feat) == 2:
            info[feat[0]] = feat[1]
        else:
            info[feat[0]] = ""

    vep_out = info["CSQ"].split(",")
    effects = [effect.split("|") for effect in vep_out]

    # Add annotation (from GFF)
    for effect in effects:
        gene = effect[4]
        if gene != "":
            gannot = annot_dict[gene]
        else:
            gannot = ""
        effect.append(gannot)

    parsed_variant = [parsed_vcf + effect for effect in effects]

    return(parsed_variant)

def parse_vep(vep_file, gff_ref, sample_names, outfld):

    gff_ref = pb.BedTool(gff_ref)
    vep = pb.BedTool(vep_file)

    # Create dict for annotation (from GFF)
    gene_types = ["gene", "pseudogene"]
    gff_gene = gff_ref.filter(lambda x: x[2] in gene_types)
    annot = {}
    for entry in gff_gene:
        ga = getAnnot(entry)
        annot[ga[0]] = ga[1]

    #sample_names = ['sample_A', 'sample_B']
    # Create DF
    col_pfx = ["RefCount_", "AltCount_", "RefRatio_", "depth_"]
    sample_cols = [[cp+n for cp in col_pfx] for n in sample_names]
    unnested = [item for sublist in sample_cols for item in sublist]
    colnames1 = ["Chrom", "Pos", "Ref", "Alt"]
    colnames2 = ["Allele", "Consequence", "IMPACT",
                 "SYMBOL", "Gene", "Feature_type",
                "Feature", "BIOTYPE", "EXON", "INTRON",
                "HGVSc", "HGVSp", "cDNA_position",
                "CDS_position", "Protein_position",
                "Amino_acids", "Codons", "Existing_variation",
                "DISTANCE", "STRAND", "FLAGS",
                "SYMBOL_SOURCE", "HGNC_ID", "SOURCE",
                "",
                "Annot"]
    colnames = colnames1+unnested+colnames2


    parsed = [parse_variant(var, annot) for var in vep]
    flat = list(chain.from_iterable(parsed))
    var_df = pd.DataFrame.from_records(flat, columns=colnames)
    outname = Path(outfld).joinpath(Path(vep_file).stem+'_parsed.tsv')
    var_df.to_csv(outname, sep = '\t')

def input_parser(input_file):
    with open(input_file, 'r+') as infile:
        lines = [l.strip() for l in infile if len(l.strip()) > 0 and not l.startswith('#')]
        pipe_dict = defaultdict(list)
        outdict = {}

        for line in lines:
            if line.startswith('>'):
                entry = line.replace('>', '')
            else:
                pipe_dict[entry].append(line)

        list_keys = ['Bams', 'Sample_names']

        for k, v in pipe_dict.items():
            if k in list_keys:
                outdict[k] = v
            else:
                outdict[k] = v[0]

        return(outdict)

#### MAIN FUNCTION ####

def pipe_main(input_file):

    print(input_file)

    start = time.time()

    ## Greet User
    print((
        '\n'
        '###################################################################\n'
        '####      Running Variant-Calling Pipeline by Lucas M.T.       ####\n'
        '###################################################################\n'
        '\n'
    ))

    #os.chdir('/mnt/Disc4T/Projects/Miniprojects/ChIP_Seq_Pipeline_Executable/')
    #input_file = './pipeline_parameters.txt'
    input_d = input_parser(input_file)
    os.makedirs(input_d['Out_folder'], exist_ok=True)
    #print(input_d)

    ## Set steps to perform
    readgroups = True if input_d['Run_readgroups'] == 'yes' else False
    markdupl = True if input_d['Run_markdupl'] == 'yes' else False
    baserecal = True if input_d['Run_baserecal'] == 'yes' else False
    appbqsr = True if input_d['Run_appbqsr'] == 'yes' else False
    mergeb = True if input_d['Run_mergeb'] == 'yes' else False
    callvars = True if input_d['Run_callvars'] == 'yes' else False
    callvep = True if input_d['Run_callvep'] == 'yes' else False
    parsevep = True if input_d['Run_parsevep'] == 'yes' else False


    steps = [k for k, v in input_d.items() if k.startswith('Run') and v == 'yes']
    print('Running steps:\n')
    for step in steps:
        print(step)

    ## Create output folders
    os.makedirs(input_d['Out_folder'], exist_ok=True)
    bamsdir = Path(input_d['Out_folder']).joinpath('Bams')
    os.makedirs(bamsdir, exist_ok=True)

    ## Add or replace readsgroups
    if readgroups:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                  Running AddOrReplaceReadGroups\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        for b in input_d['Bams']:
            bpath = Path(input_d['In_fld']).joinpath(b)
            AddOrReplaceReadGroups(input_d['GATK_path'], bpath, bamsdir)

    ## Mark duplicates
    if markdupl:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                    Running MarkDuplicates\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        for b in input_d['Bams']:
            bpath = bamsdir.joinpath(Path(b).stem+'_withRG.bam')
            mark_duplicates(input_d['GATK_path'], bpath, bamsdir)

    ## Base recalibration
    if baserecal:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                    Running BaseRecalibrator\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        for b in input_d['Bams']:
            bpath_1 = bamsdir.joinpath(Path(b).stem+'_withRG.bam')
            bpath_2 = str(bpath_1).replace('.bam', '_markedDuplicates.bam')
            base_recalibration(input_d['GATK_path'],
                               bpath_2, bamsdir,
                               input_d['Ref_bed'],
                               input_d['Ref_fa'])

    ## Apply Base Quality Score Recalibration
    if appbqsr:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                    Running ApplyBQSR\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        for b in input_d['Bams']:
            bpath_1 = bamsdir.joinpath(Path(b).stem+'_withRG.bam')
            bpath_2 = str(bpath_1).replace('.bam', '_markedDuplicates.bam')
            applyBQSR(input_d['GATK_path'],
                      bpath_2, bamsdir,
                      input_d['Ref_fa'])

    ## Merge Bams
    if mergeb:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                         Merging Bams\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        inbams = []
        for b in input_d['Bams']:
            bpath_1 = bamsdir.joinpath(Path(b).stem+'_withRG.bam')
            bpath_2 = str(bpath_1).replace('.bam', '_markedDuplicates_BQSR.bam')
            inbams.append(bpath_2)

        outbam = bamsdir.joinpath('all_samples_merged')
        mergeBams(input_d['GATK_path'], *inbams,
                  outfile = outbam)
        samtools_index(input_d['Sam_path'], str(outbam)+'.bam')

    ## Merge Bams
    if callvars:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                         Calling Variants\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        mergedbam = Path(input_d['Out_folder']).joinpath('Bams/all_samples_merged.bam')
        call_variants(input_d['GATK_path'],
                      input_d['Ref_fa'],
                      mergedbam,
                      Path(input_d['Out_folder']).joinpath('all_samples_variants.vcf'))

    ## Annotate variants with VEP
    if callvep:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                Annotating Variants with VEP\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        vcf = str(Path(input_d['Out_folder']).joinpath('all_samples_variants.vcf'))
        call_VEP(input_d['Vep_path'], vcf,
                 input_d['Ref_gff'],
                 input_d['Ref_fa'])

    ## Annotate variants with VEP
    if parsevep:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                    Parsing VEP annotations\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        ## Parse VEP variants
        vep_name = 'all_samples_variants_VEPannotated.txt'
        vep_file = str(Path(input_d['Out_folder']).joinpath(vep_name))
        parse_vep(vep_file, input_d['Ref_gff'],
                  input_d['Sample_names'], input_d['Out_folder'])



    ## Time-it and Finish

    end = time.time()
    out_fld = input_d['Out_folder']
    str_time = time.strftime("%H:%M:%S", time.gmtime(end - start))

    print((
        '\n'
        '###################################################################\n'
        '##                                                               ##\n'
        '##                  Finished!                                    ##\n'
        f'##                  Elapsed time: {str_time}                       ##\n'
        '##                                                               ##\n'
        '###################################################################\n'
    ))

#### Parse Arguments and Run ####


def run():
    '''
    Parse command line args and run 'pipe_main(args)'.
    '''
    program_description = (
        'Pipeline for ChIP-Seq Data analysis:'
        'Takes the \'pipeline_parameters.txt\' file as input '
        'and performs all the steps necessary for the analysis.'
    )

    parser = argparse.ArgumentParser(description=program_description)

    # Required Arguments
    hline = (
        'Parameters file: all the parameters for the analysis '
        'must be set in this text file. Follow '
        'the instructions in it!'
    )
    parser.add_argument('-pf', type = os.path.abspath, dest = 'params',
                        metavar = 'params_file', required = True,
                        help=hline)

    args = parser.parse_args()

    pipe_main(
        input_file = args.params
    )

if __name__ == "__main__":
    run()


