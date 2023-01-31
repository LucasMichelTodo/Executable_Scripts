import os
import subprocess as sp
import multiprocessing as mp
import argparse

## Functions

def subset_fastq(fastq, nreads, np):
    ext = '.'+'.'.join(fastq.split('.')[1:])
    out_fq = fastq.replace(ext, f'_subset{nreads}.fq')

    cmd = ['seqtk', 'sample', '-s123', fastq, str(nreads)]
    print(f'Subseting {fastq} to {out_fq} ({str(nreads)} reads) ...')
    subset = sp.run(cmd, stdout=sp.PIPE).stdout.decode('utf-8')
    with open(out_fq, 'w+') as outfile:
        outfile.write(subset)

def subset_fastq_paired_parallel(fastq1, fastq2, nreads, np):
    pool = mp.Pool(np)
    pool.starmap(subset_fastq, [[fastq1, nreads, np], [fastq2, nreads, np]])

## Parse Arguments and run program

wd = os.path.dirname(os.path.realpath(__file__))

def run():
    '''
    Parse command line args and run 'subset_fastq_paired_parallel(args)'.
    '''
    program_description = ('\'Fastq\' subset maker: '
                           'It takes as input the reads1 and reads2 files of a '
                           'paired-end -Seq file, a number of reads and a number  '
                           'of CPU processors and it returns randomly subseted \'.fq\' files '
                           'with the desired number of reads only.'
                           )

    parser = argparse.ArgumentParser(description=program_description)

    # Required Arguments

    hline = 'Reads file in \'fastq\' format 1'
    parser.add_argument('-r1', type = str, dest = 'reads_1',
                        metavar = 'reads_1', required = True,
                        help=hline)

    hline = 'Reads file in \'fastq\' format 2'
    parser.add_argument('-r2', type = str, dest = 'reads_2',
                        metavar = 'reads_2', required = True,
                        help=hline)

    hline = 'Number of total reads for the subseted files'
    parser.add_argument('-n', type = int, dest = 'nreads',
                        metavar = 'num_reads', required = True,
                        help=hline)

    # Optional Arguments
    parser.add_argument('-np', type = int, dest = 'threads',
                        metavar = 'num_process', default = 1,
                        help = 'Number of processors to use.')

    args = parser.parse_args()

    subset_fastq_paired_parallel(
        fastq1 = args.reads_1,
        fastq2 = args.reads_2,
        nreads = args.nreads,
        np = args.threads
    )


if __name__ == "__main__":
    run()


########################################
# wd = '/mnt/Disc4T/Projects/Alba/Mostres_13_09_22/Raw_Data/'
# os.chdir(wd)

# fastq1s = [f for f in os.listdir()]
# fastq2s = [
#     'NCV26_lib_08246AAC_ATCACG_R2_001.fastq.gz',
#     'NCV27_lib_08247AAC_CGATGT_R2_001.fastq.gz',
#     'NCV28_lib_08248AAC_TTAGGC_R2_001.fastq.gz',
#     'NCV29_lib_08249AAC_TGACCA_R2_001.fastq.gz',
# ]

# for f1, f2 in zip(fastq1s, fastq2s):
#     subset_fastq_paired_parallel(f1, f2, 1000)
