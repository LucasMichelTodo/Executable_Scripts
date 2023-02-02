import os
import sys
import subprocess as sp
import pandas as pd
import re
from datetime import datetime
from collections import defaultdict
import argparse
import time

#### FUNCTIONS ####

## Pre-Processing

def get_extension(filename):
    '''
    Get extension from a filename. It can have up to 2 extensions,
    a filetype extension and optionally a compressing extension (e.g., '.gz')
    '''
    p = pathlib.Path(filename)
    s = p.suffixes[-2:]
    return(''.join(s))

def call_BBDUK(program_path, in1, in2, ref, outfolder='./', params=''):

    if not outfolder.endswith('/'): outfolder = outfolder+'/'
    extension = get_extension(in1)

    out1 = outfolder + in1.rsplit("/", 1)[1].replace(extension, "_clean.fq")
    out2 = outfolder + in2.rsplit("/", 1)[1].replace(extension, "_clean.fq")
    outm = outfolder + in1.rsplit("/", 1)[1].replace(extension, "_badreads.fq")

    cmd = [
        program_path,
        'in='+in1, 'in2='+in2,
        'out='+out1, 'out2='+out2, 'outm='+outm,
        'ref='+ref,
    ]

    params = params.split(' ')
    cmd = cmd+params
    sp.run(cmd)

def call_fastq_screen(program_path, fastq_file_list, threads = '1', outdir = './'):
    ## input a list of fastq files to be screened
    cmd = [program_path, '--threads', threads, '--outdir', outdir] + fastq_file_list
    sp.run(cmd)

def call_fastqc(program_path, fastq_file_list, threads = '1', outdir = './'):
    ## input a list of fastq files to be screened
    cmd = [program_path, '--threads', threads, '--outdir', outdir] + fastq_file_list
    sp.run(cmd)

def call_multiqc(program_path, fld):
    cmd = [program_path, fld, '-f', '-o', f'{fld}MultiQC/']
    sp.run(cmd)

## Alignment

def call_Bowtie2(program_path, in1, in2, out, threads, params):
    cmd = [program_path, '-1', in1, '-2', in2, '-S', out, '-p', threads]
    params = params.split(' ')
    rep_name = out.rsplit('.', 1)[0]+'_bwt2_align_report.txt'

    cmd = cmd+params
    print(cmd)
    with open(rep_name, 'w+') as outfile:
        outfile.write(f'{out}\n\n'+' '.join(cmd)+'\n')

    result = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
    str_result = result.stdout.decode('utf-8')

    with open(rep_name, 'a+') as outfile:
        outfile.write(str_result)

def parse_agregated_report(report, start_line, out_name):
    with open(report, 'r+') as infile:

        infile = open(report, 'r+')
        lines = infile.readlines()
        lines = [line.rstrip() for line in lines]
        out = {}

        for idx, line in enumerate(lines):

            if line.startswith(start_line):

                name = line.split('/')[-1].strip()
                out[name] = {}

                out[name]['call'] = lines[idx+2]
                out[name]['nreads'] = lines[idx+3].split()[0]

                paired = re.search('\((.*)\)', lines[idx+4])
                out[name]['paired'] = paired.group(1)
                conc_0 = re.search('\((.*)\)', lines[idx+5])
                out[name]['conc_0'] = conc_0.group(1)
                conc_1 = re.search('\((.*)\)', lines[idx+6])
                out[name]['conc_1'] = conc_1.group(1)
                conc_more = re.search('\((.*)\)', lines[idx+7])
                out[name]['conc_more'] = conc_more.group(1)

                conc_0_disc1 = re.search('\((.*)\)', lines[idx+10])
                out[name]['conc_0_disc1'] = conc_0_disc1.group(1)

                conc_0_disc0_al0 = re.search('\((.*)\)', lines[idx+14])
                out[name]['conc_0_disc0_al0'] = conc_0_disc0_al0.group(1)
                conc_0_disc0_al1 = re.search('\((.*)\)', lines[idx+15])
                out[name]['conc_0_disc0_al1'] = conc_0_disc0_al1.group(1)
                conc_0_disc0_almore = re.search('\((.*)\)', lines[idx+16])
                out[name]['conc_0_disc0_almore'] = conc_0_disc0_almore.group(1)

                out[name]['overall_al'] = lines[idx+17].replace(' overall alignment rate', '')

        out_df = pd.DataFrame.from_dict(out, orient='index')
        # Reorder cols
        cols = out_df.columns.tolist()
        cols = [cols[1]]+[cols[-1]]+cols[2:]+[cols[0]]
        out_df = out_df[cols]

        new_cols = [
            'Total Reads',
            'Overall Alignment Rate',
            'Paired Reads',
            'Aligned Concordantly 0 times',
            'Aligned Concordantly 1 time',
            'Aligned Concordantly >1 times',
            'Aligned Concordantly 0 times, Discordantly 1',
            'Non-concordant or discordant, Mate Aligned 0 times',
            'Non-concordant or discordant, Mate Aligned 1 time',
            'Non-concordant or discordant, Mate Aligned >1 times',
        ]

        re_cols = dict(zip(cols, new_cols))
        out_df.rename(columns=re_cols, inplace=True)
        out_df.to_csv(out_name, sep='\t')

def common_start(sa, sb):
    """ returns the longest common substring from the beginning of sa and sb """
    def _iter():
        for a, b in zip(sa, sb):
            if a == b:
                yield a
            else:
                return

    return ''.join(_iter())

def bowtie2_align_pipe(bowtie2path, sampath, read1s_list, read2s_list, names, threads, params, outpath = './'):
    """
    Take a sorted list of read1, read 2 files and their desired filenames
    (all lists must be in the same sort order), a params string with all wanted
    bowtie2 parameters and optionally a path for the outputs. Runs the alignments
    and generates an automated report.
    """
    startTime = datetime.now()
    os.makedirs(outpath, exist_ok=True)

    for r1, r2, name in zip(read1s_list, read2s_list, names):
        out = outpath+name+'.sam'
        print('\n'+f'Aligning: {out}'+'\n')
        call_Bowtie2(bowtie2path, r1, r2, out, threads, params)
        from_sam_to_bam(sampath, out, threads)

    finishTime = datetime.now()
    time_lines = (
        f'\nProcess started at {startTime}\n'
        f'Process ended at {finishTime}\n'
        f'Process started at {finishTime - startTime}'
    )

    ## Write report
    reports = [f for f in os.listdir(outpath) if f.endswith('bwt2_align_report.txt')]
    with open(outpath+'full_report.txt', 'w+') as outfile:
        for fname in reports:
            with open(outpath+fname) as infile:
                outfile.write(infile.read())
        outfile.write(time_lines)

    ## Generate table report
    parse_agregated_report(outpath+'full_report.txt', outpath, outpath+'full_report_table.tsv')


## Samtools
def from_sam_to_bam(program_path, samfile, threads):

    name = samfile.rsplit(".", 1)[0]
    cmd = "{} view -@ {} -bS {} > {}" .format(program_path, threads, samfile, name+".bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ### Erase SAM after creating BAM
    cmd = "rm {}" .format(samfile)
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "{} sort -@ {} {} > {}" .format(program_path, threads, name+".bam", name+"_sort.bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ### Erase bam after creating sortedBAM
    cmd = "rm {}" .format(name+".bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "{} index -@ {} {} > {}" .format(program_path, threads, name+"_sort.bam", name+"_sort.bam.bai")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ## Filter only >=q5 reads
    cmd = "{} view -@ {} -b -q 5 {} > {}" .format(program_path, threads, name+"_sort.bam", name+"_q5_sort.bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "{} index -@ {} {} > {}" .format(program_path, threads, name+"_q5_sort.bam", name+"_q5_sort.bam.bai")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

## Remove Duplicates

def remove_duplicates(program_path, outdir, bam):

    name = bam.rsplit("/", 1)[1]
    o = outdir+name.replace(".bam", "_noDup.bam")
    m = outdir+name.replace(".bam", "_metrics.txt")

    cmd = (f"{program_path} MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(bam, o, m)

    sp.call(cmd, shell=True)

def samtools_index(program_path, bam, np=8):
    cmd = [program_path, 'index', '-@', str(np), bam]
    print(' '.join(cmd))
    sp.run(cmd)

## Get coverage (Tracks)

def get_RPKMs(program_path, bam, outfld, threads, params):

    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam.rsplit('/', 1)[1]
    outname = outfld+name.replace('.bam', f'.bdg')

    cmd = (f'{program_path} -b {bam} -o {outname} -p {threads} {params}')
    sp.call(cmd, shell=True)


def get_RPKMs_normInput(program_path, bam_IP, bam_in, outfld, threads, params):

    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam_IP.rsplit('/', 1)[1]
    outname = outfld+name.replace('.bam', '_normInput.bdg')

    cmd = (
        f'{program_path} -b1 {bam_IP} -b2 {bam_in} '
        f'-o {outname} -p {threads} {params}'
    )
    sp.call(cmd, shell=True)

def params_to_name(params_dict):
    """
    Take a dicttionary with params and convert it to a string to add to name files.
    """
    suffix = ''
    for k,v in params_dict.items():
        k_str = k.replace('-', '').strip()
        v_str = v.strip()
        suffix += f'_{k_str}{v_str}'
    return(suffix)

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

        list_keys = ['Read1s', 'Read2s', 'Sample_names']

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
        '####          Running RNA-Seq Pipeline by Lucas M.T.           ####\n'
        '###################################################################\n'
        '\n'
    ))

    input_d = input_parser(input_file)
    os.makedirs(input_d['Out_folder'], exist_ok=True)
    input_d = input_parser(input_file)
    if not input_d['Reads_folder'].endswith('/'):
        input_d['Reads_folder'] = input_d['Reads_folder']+'/'
    #print(input_d)

    ## Set steps to perform
    clean = True if input_d['Run_BBDUK'] == 'yes' else False
    fastqc = True if input_d['Run_FastQC'] == 'yes' else False
    fastscreen = True if input_d['Run_Fastq_Screen'] == 'yes' else False
    align = True if input_d['Run_Bowtie2'] == 'yes' else False
    deduplicate = True if input_d['Run_GATK'] == 'yes' else False
    tracks = True if input_d['Run_BamCoverage'] == 'yes' else False

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

        print('Cleaning Reads... ')
        for f in input_d['Read1s']: print(f)
        print('-----')
        for f in input_d['Read2s']: print(f)
        print('-----')

        outpath = input_d['Out_folder']+'/Clean_Reads/'
        os.makedirs(outpath, exist_ok=True)

        read1s = [f for f in input_d['Read1s']]
        read2s = [f for f in input_d['Read2s']]

        for pair in zip(read1s, read2s):
            print(pair)
            in1 = input_d['Reads_folder']+pair[0]
            in2 = input_d['Reads_folder']+pair[1]
            call_BBDUK(
                input_d['BBDUK_path'], in1, in2,
                input_d['BBDUK_ref'], outpath, input_d['BBDUK_params']
            )

    ## Call FastQC
    if fastqc:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                Calling FastQC on clean reads\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        cl_dir = input_d['Out_folder']+'/Clean_Reads/'
        all_reads = [cl_dir+f for f in os.listdir(cl_dir) if f.endswith('clean.fq')]
        outdir = input_d['Out_folder']+'/FastQC_Clean/'
        os.makedirs(outdir, exist_ok=True)
        call_fastqc(
            input_d['FastQC_path'], all_reads,
            threads = input_d['Threads'], outdir = outdir
        )
        call_multiqc(input_d['MultiQC_path'], outdir)

    ## Call fastq_screen on clean reads
    if fastscreen:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '             Calling Fastq_Screen on clean reads\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        cl_dir = input_d['Out_folder']+'/Clean_Reads/'
        all_reads = [cl_dir+f for f in os.listdir(cl_dir) if f.endswith('clean.fq')]
        outdir = input_d['Out_folder']+'/Fastq_Screens_Clean/'
        os.makedirs(outdir, exist_ok=True)
        call_fastq_screen(
            input_d['Fastq_Screen_path'], all_reads,
            threads = input_d['Threads'], outdir=outdir
        )
        call_multiqc(input_d['MultiQC_path'], outdir)

    ## Algin Clean Reads
    if align:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                  Aligning clean reads\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        indir = input_d['Out_folder']+'/Clean_Reads/'
        outdir = input_d['Out_folder']+'/Alignments/'
        os.makedirs(outdir, exist_ok=True)

        all_reads1 = input_d['Read1s']
        all_reads2 = input_d['Read2s']
        all_names = input_d['Sample_names']
        reads1 = [indir+f.replace(get_extension(f), '_clean.fq') for f in all_reads1]
        reads2 = [indir+f.replace(get_extension(f), '_clean.fq') for f in all_reads2]

        bowtie2_align_pipe(
            input_d['Bowtie2_path'], input_d['Samtools_path'],
            reads1, reads2, all_names,
            input_d['Threads'], input_d['Bowtie2_params'], outdir
        )

    ## Remove Duplicates (Optional Step)
    if deduplicate:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                     Removing Duplicates\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        indir = input_d['Out_folder']+'/Alignments/'
        all_names = input_d['Sample_names']
        fnames = [n+"_q5_sort.bam" for n in all_names]
        outdir = indir+'DeDuplicated/'
        os.makedirs(outdir, exist_ok=True)

        for fname in fnames:
            remove_duplicates(input_d['GATK_path'], outdir, indir+fname)

        for fname in fnames:
            infile = outdir+fname.replace('.bam', '_noDup.bam')
            samtools_index(input_d['Samtools_path'], infile, input_d['Threads'])

    ## Make tracks
    if tracks:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                      Creating Tracks\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        indir = input_d['Out_folder']+'/Alignments/'
        indir = indir+'DeDuplicated/' if deduplicate else indir
        suffix = '_q5_sort_noDup.bam' if deduplicate else '_q5_sort.bam'
        all_samples = [name+suffix for name in input_d['Sample_names']]
        outdir = input_d['Out_folder']+'/Tracks/'
        os.makedirs(outdir, exist_ok=True)

        for s in all_samples:
            get_RPKMs(
                input_d['BamCoverage_path'], indir+s,
                outdir, input_d['Threads'],
                input_d['BamCoverage_params']
            )
            
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

