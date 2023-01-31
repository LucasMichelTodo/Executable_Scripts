#### Imports ####

import subprocess as sp
import os

gatk = '/home/lucas/Programs/gatk-4.1.9.0/gatk'
wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

### FUNTIONS ####

def index_feature_file(feat_file):
    cmd = '{} IndexFeatureFile -I {}' .format(gatk, feat_file)
    sp.call(cmd, shell = True)

def AddOrReplaceReadGroups(bam):
    output = bam.replace('.bam', '_withRG.bam')
    name = bam.rsplit('.')[0]
    cmd = ("java -jar /home/lucas/Programs/picard.jar "
           "AddOrReplaceReadGroups "
           "-INPUT {} "
           "-OUTPUT {} "
           "-RGID group_{} "
           "-RGLB lib_{} "
           "-RGPL illumina "
           "-RGPU unit1 "
           "-RGSM {}_sample") .format(bam, output, name, name, name)
    sp.call(cmd, shell=True)

def mark_duplicates(bam):
    outfile = bam.replace('.bam', '_markedDuplicates.bam')
    mtrcsfile = bam.replace('.bam', '_metrics.txt')
    args = [gatk, bam, outfile, mtrcsfile]
    cmd = '{} MarkDuplicates -I {} -O {} -M {}' .format(*args)
    sp.call(cmd, shell=True)

    cmd = 'samtools sort {} -o {}' .format(outfile, outfile)
    sp.call(cmd, shell=True)

def base_recalibration(bam):

    outfile = bam.replace('.bam', '_baserecal_table.table')
    known = './empty.bed'
    ref = './ref.fasta'
    args = [gatk, bam, ref, known, outfile]

    cmd = ('{} BaseRecalibrator '
           '-I {} -R {} '
           '--known-sites {} '
           '-O {}') .format(*args)

    sp.call(cmd, shell=True)

def applyBQSR(bam):

    outfile = bam.replace('.bam', '_BQSR.bam')
    recal_table = bam.replace('.bam', '_baserecal_table.table')
    ref = './ref.fasta'
    args = [gatk, bam, ref, recal_table, outfile]

    cmd = ('{} ApplyBQSR '
           '-I {} -R {} '
           '--bqsr-recal-file {} '
           '-O {}') .format(*args)

    sp.call(cmd, shell=True)

def mergeBams(*bams, out):
    nbams = len(bams)
    inputs = '-I {} '*nbams
    cmd = 'java -jar /home/lucas/Programs/picard.jar ' \
        'MergeSamFiles ' + \
        inputs .format(*bams) + \
        '-O {}.bam' .format(out)
    sp.call(cmd, shell=True)

def call_variants(bam):
    outfile = bam.replace('.bam', '_variants.vcf')
    ref = './ref.fasta'
    args = [gatk, ref, bam, outfile]

    cmd = ('{} --java-options "-Xmx4g" HaplotypeCaller '
           '-R {} -I {} -O {} -ploidy 1') .format(*args)

    sp.call(cmd, shell=True)

def call_VEP(vcf, gff, fasta):

    out = vcf.replace('.vcf', '_VEPannotated.txt')
    args = [vcf, out, gff, fasta]

    cmd = ("/home/lucas/Programs/ensembl-vep/vep "
           "-i {} "
           "-o {} "
           "--gff {} "
           "--fasta {} "
           "--force_overwrite "
           "--vcf") .format(*args)

    sp.call(cmd, shell=True)

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


def pipe_main(input_file):

    start = time.time()

    ## Greet User
    print((
        '\n'
        '###################################################################\n'
        '####      Running Variant-Calling Pipeline by Lucas M.T.       ####\n'
        '###################################################################\n'
        '\n'
    ))

    input_d = input_parser(input_file)
    os.makedirs(input_d['Out_folder'], exist_ok=True)
    input_d = input_parser(input_file)
    #print(input_d)

    ## Set steps to perform
    clean = True if input_d['Run_BBDUK'] == 'yes' else False
    fastqc = True if input_d['Run_FastQC'] == 'yes' else False
    fastscreen = True if input_d['Run_Fastq_Screen'] == 'yes' else False
    align = True if input_d['Run_Bowtie2'] == 'yes' else False
    deduplicate = True if input_d['Run_GATK'] == 'yes' else False
    rawtracks = True if input_d['Run_BamCoverage'] == 'yes' else False
    normtracks = True if input_d['Run_BamCompare'] == 'yes' else False

    steps = [k for k, v in input_d.items() if k.startswith('Run') and v == 'yes']
    print('Running steps:\n')
    for step in steps:
        print(step)

    ## Clean Reads
    if clean:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                       Cleaning Reads\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

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
        f'##                  Results in: {out_fld}                   ##\n'
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
    parser.add_argument('-pf', type = str, dest = 'params',
                        metavar = 'params_file', required = True,
                        help=hline)

    args = parser.parse_args()

    pipe_main(
        input_file = args.params
    )

if __name__ == "__main__":
    run()



####################################################

import os
import subprocess as sp

wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

gatk = '/home/lucas/Programs/gatk-4.1.9.0/gatk'
indir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Bams/'

os.listdir(indir)

bams = ['1.2B_in_sort_q5.bam',
        '10G_in_sort_q5.bam',
        'A7K9_in_sort_q5.bam',
        'E5K9_in_sort_q5.bam',
        'B11_in_sort_q5.bam'
        #'NF54_in_renamed_q5_sort.bam',
        ]

index_feature_file('./empty.bed')

for bam in bams:
    bam = indir+bam

    AddOrReplaceReadGroups(bam)
    bam = bam.replace('.bam', '_withRG.bam')

    mark_duplicates(bam)
    bam = bam.replace('.bam', '_markedDuplicates.bam')

    base_recalibration(bam)
    applyBQSR(bam)
    bam = bam.replace('.bam', '_BQSR.bam')


bamlist = [f for f in os.listdir(indir) if f.endswith('_withRG_markedDuplicates_BQSR.bam')]

os.chdir(indir)
mergeBams(*bamlist, out = 'merged_12B_10G_A7_E5_B11')

bam = 'merged_12B_10G_A7_E5_B11.bam'
sp.call('samtools index {}' .format(bam), shell=True)
call_variants(bam)

os.chdir(wd)
vcf = './merged_12B_10G_A7_E5_B11_variants.vcf'
gff = './PlDB-52_Pfalciparum3D7_vep_changetypes.gff.gz'
fasta = './ref.fasta'
call_VEP(vcf, gff, fasta)
