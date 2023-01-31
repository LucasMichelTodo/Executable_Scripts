import time
import os
import pybedtools as pb
import numpy as np
import subprocess as sp
from tqdm import tqdm
from sklearn.mixture import GaussianMixture
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import matplotlib.pyplot as plt
import pandas as pd
from itertools import repeat
import argparse


## Functions

def create_window_ref(window_size, genome_file, step):
    bed = pb.BedTool()
    win_ref = bed.window_maker(w=window_size, g = genome_file, s = step)
    return(win_ref)

def plot_gaussian_mixture(cov_vector, gm, plotname):

    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    logprob = gm.score_samples(x.reshape(-1, 1))
    pdf = np.exp(logprob)
    plt.plot(x, pdf, 'k', linewidth=2)
    plt.savefig(plotname)
    plt.clf()

def plot_components(cov_vector, gm, plotname):
    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    logprob = gm.predict_proba(x.reshape(-1, 1))
    cdf1 = logprob[:,0]
    cdf2 = logprob[:,1]
    plt.plot(x, cdf1, 'k', linewidth=2, color='blue', label="Component 1")
    plt.plot(x, cdf2, 'k', linewidth=2, color='green', label="Component 2")
    plt.legend(loc='upper left')
    plt.savefig(plotname)
    plt.clf()


def get_probs(coverage, gm):
    return(gm.predict_proba(coverage.reshape(-1, 1)))

def get_differential_peaks(
        peakfile1, peakfile2,
        covfile1, covfile2,
        prefix1, prefix2,
        genome_file,
        window_size, stepsize,
        minprobdif,
        mergedist,
        minlen,
        outfld,
        num_p):

    start = time.time()

    ## Greet User
    print((
        '\n'
        '###################################################################\n'
        '###   Running Custom differential Peak-Calling by Lucas M.T.    ###\n'
        '###################################################################\n'
        '\n'
    ))

    ## Create folders for output
    os.makedirs(outfld, exist_ok = True)
    if not outfld.endswith('/'): outfld = outfld+'/'
    os.makedirs(outfld+'Data/', exist_ok = True)
    os.makedirs(outfld+'Data/Common_Coverage/', exist_ok = True)
    os.makedirs(outfld+'Data/Common_Peaks/', exist_ok = True)
    os.makedirs(outfld+'Data/Plots/', exist_ok = True)
    os.makedirs(outfld+'Data/Window_Coverage/', exist_ok = True)
    os.makedirs(outfld+'Data/PreFilter_Difpeaks/', exist_ok = True)


    ## Generate GaussianMixture distribution with 2 components for coverage in peaks
    ## Peaks coverage
    print('Fitting distribution...')

    pf1 = pb.BedTool(peakfile1)
    pf2 = pb.BedTool(peakfile2)
    cf1 = pb.BedTool(covfile1)
    cf2 = pb.BedTool(covfile2)

    cov_peaks1 = np.array([float(f.fields[10]) for f in pf1.map(cf1, c = 4, o='mean')])
    cov_peaks2 = np.array([float(f.fields[10]) for f in pf2.map(cf2, c = 4, o='mean')])

    ## All sklearn estimators use as input a 2D array,
    ## with samples as rows and features as columns.
    ## Our data-set consists in one sample per bed feature
    ## with only one feature (coverage), so we have to reshape it
    ## into a 2D array with just one column.
    ## -1 in reshape just means whatever it takes to make it work!
    ## So (-1, 1) means: any number of rows and 1 column

    gm1 = GaussianMixture(n_components=2).fit(cov_peaks1.reshape(-1, 1))
    gm2 = GaussianMixture(n_components=2).fit(cov_peaks2.reshape(-1, 1))

    suffix = '_peak_coverage_fit.png'
    plot_gaussian_mixture(cov_peaks1, gm1, f'{outfld}Data/Plots/{prefix1}{suffix}')
    plot_gaussian_mixture(cov_peaks2, gm2, f'{outfld}Data/Plots/{prefix2}{suffix}')
    plot_components(cov_peaks1, gm1, f'{outfld}Data/Plots/{prefix1}_components_cdf.png')
    plot_components(cov_peaks2, gm2, f'{outfld}Data/Plots/{prefix2}_components_cdf.png')

    ## Create window reference
    print('Creating windowed coverage...')
    suffix = '_window_coverage.bdg'
    win_ref = create_window_ref(window_size, genome_file, stepsize)
    win_cov1 = win_ref.map(cf1, c=4, o='mean').saveas(f'{outfld}Data/Window_Coverage/{prefix1}{suffix}')
    win_cov2 = win_ref.map(cf2, c=4, o='mean').saveas(f'{outfld}Data/Window_Coverage/{prefix2}{suffix}')

    ## Create join coverage bdg
    print('Creating join coverage bed...')
    outfile = outfld+f'Data/Common_Coverage/{prefix1}_{prefix2}_common_coverage.bdg'
    cmd = ['bedtools', 'unionbedg', '-i',
           f'{outfld}Data/Window_Coverage/{prefix1}{suffix}',
           f'{outfld}Data/Window_Coverage/{prefix2}{suffix}',
           '>', outfile]
    sp.call(' '.join(cmd), shell = True)
    union_bed = pb.BedTool(outfile)

    ## Create common peaks bed
    print('Creating common peaks bed...')

    files_to_cross = []
    str_beds = peakfile1 + ' ' + peakfile2
    outfile = f'{outfld}Data/Common_Peaks/{prefix1}_{prefix2}_common_peaks.bed'
    cmd = 'awk \'{print}\' '+f'{str_beds} > {outfile}'
    sp.call(cmd, shell = True)
    common_peaks = pb.BedTool(outfile).sort()

    ## Subset union_bed to common_peaks
    common_cov_common_peaks = union_bed.intersect(common_peaks)

    ## Call differential peaks
    print('Calling differential peaks (this might take a while)...')

    c1s = [float(x.fields[3]) for x in common_cov_common_peaks]
    c2s = [float(x.fields[4]) for x in common_cov_common_peaks]

    c1s = np.array(c1s)
    c2s = np.array(c2s)

    prob1s_pre = get_probs(c1s, gm1)
    prob2s_pre = get_probs(c2s, gm2)
    prob1s = prob1s_pre[:,0]
    prob2s = prob2s_pre[:,0]

    probdifs = prob1s - prob2s
    peaks1over2 = probdifs > minprobdif
    peaks2over1 = -probdifs > minprobdif
    chroms = [x.chrom for x in common_cov_common_peaks]
    starts = [x.start for x in common_cov_common_peaks]
    stops = [x.stop for x in common_cov_common_peaks]

    pd.DataFrame(common_cov_common_peaks)
    df = pd.read_table(
        common_cov_common_peaks.fn,
        names=['#chrom', 'start', 'stop', 'cov1', 'cov2']
               )
    rawout1 = (f'{outfld}Data/PreFilter_Difpeaks/'
               f'difpeaks_{prefix1}over{prefix2}'
               f'_w{window_size}_s{stepsize}_pd{minprobdif}.bed')

    rawout2 = (f'{outfld}Data/PreFilter_Difpeaks/'
               f'difpeaks_{prefix2}over{prefix1}'
               f'_w{window_size}_s{stepsize}_pd{minprobdif}.bed')

    df[peaks1over2].to_csv(rawout1, sep='\t', index=False)
    df[peaks2over1].to_csv(rawout2, sep='\t', index=False)

    rawbed1 = pb.BedTool(rawout1)
    rawbed2 = pb.BedTool(rawout2)

    ## Merge peaks and filter by length
    out1 = (f'{outfld}difpeaks_{prefix1}_over_{prefix2}'
            f'_w{window_size}_s{stepsize}_pd{minprobdif}'
            f'_mg{mergedist}_ml{minlen}'
            '_difpeaks.bed')

    out2 = (f'{outfld}difpeaks_{prefix2}_over_{prefix1}'
            f'_w{window_size}_s{stepsize}_pd{minprobdif}'
            f'_mg{mergedist}_ml{minlen}'
            '_difpeaks.bed')

    out1_pre = rawbed1.sort().merge(d=mergedist)
    out1 = out1_pre.filter(lambda f: f.stop - f.start > minlen).saveas(out1)
    out2_pre = rawbed2.sort().merge(d=mergedist)
    out2 = out2_pre.filter(lambda f: f.stop - f.start > minlen).saveas(out2)

    end = time.time()
    print('Finished! Elapsed time:')
    print(end - start)

## Parse Arguments and run program

wd = os.path.dirname(os.path.realpath(__file__))

def run():
    '''
    Parse command line args and run 'get_differential_peaks(args)'.
    '''
    program_description = ('Differential Peaks Caller:'
                           'For each sample, takes as input a \'.narrowPeak\' file '
                           '(from a MACS2 callpeak call) '
                           'and a bed/bedgraph coverage file '
                           '(generated with DeepTools e.g.) '
                           'and returns a list of differential peaks.')

    parser = argparse.ArgumentParser(description=program_description)

    # Required Arguments

    hline = ('Peaks file in MACS2 format '
             '(.narrowPeaks) for file 1')
    parser.add_argument('-p1', type = str, dest = 'peaks_1',
                        metavar = 'peaks_1', required = True,
                        help=hline)

    hline = ('Peaks file in MACS2 format '
             '(.narrowPeaks) for file 2')
    parser.add_argument('-p2', type = str, dest = 'peaks_2',
                        metavar = 'peaks_2', required = True,
                        help=hline)

    hline = ('Coverage file in BED/BEDGRAPH format '
             'for file 1')
    parser.add_argument('-c1', type = str, dest = 'cov_1',
                        metavar = 'coverage_1', required = True,
                        help=hline)

    hline = ('Coverage file in BED/BEDGRAPH format '
             'for file 2')
    parser.add_argument('-c2', type = str, dest = 'cov_2',
                        metavar = 'coverage_2', required = True,
                        help=hline)

    hline = ('A \'.genome\' file with one gene chromosome '
             'per line and its length, separated by tabs')
    parser.add_argument('-g', type = str, dest = 'genome_file',
                        metavar = 'genome_file', required = True,
                        help = hline)

    # Optional Arguments
    parser.add_argument('-pfx1', type = str, dest = 'prefix_1',
                        metavar = 'prefix_1', default = 'sample_1',
                        help = 'Prefix for sample 1 (sample_1)')

    parser.add_argument('-pfx2', type = str, dest = 'prefix_2',
                        metavar = 'prefix_2', default = 'sample_2',
                        help = 'Prefix for sample 2 (sample_2)')

    parser.add_argument('-ws', type = int, dest = 'window_size',
                        metavar = 'window_size', default = 100,
                        help = 'Window size for coverage calculations (100)')

    hline = ('Step size for coverage calculations (100). '
             'Must be <= to window size')
    parser.add_argument('-ss', type = int, dest = 'step_size',
                        metavar = 'step_size', default = 100,
                        help = hline)

    hline = ('Minimum probability difference between '
             'sample distributions for a '
             'differential peak to be called (0.3). '
             'A bigger number means a more stringent call')
    parser.add_argument('-mpd', type = float, dest = 'min_prob_dif',
                        metavar = 'min_prob_diff', default = 0.3,
                        help = hline)

    hline = ('Differential peaks separated by less '
             'than this distance will be merged together (500)')
    parser.add_argument('-md', type = int, dest = 'merge_dist',
                        metavar = 'merge_distance', default = 500,
                        help = hline)

    hline = ('Differential peaks smaller than '
             'this length will be discarded (1000). '
             'Peaks are first joined by merge_dist and then '
             'filtered by min_len')
    parser.add_argument('-ml', type = int, dest = 'min_len',
                        metavar = 'min_length', default = 1000,
                        help = hline)

    parser.add_argument('-o', '--out_folder', type = str, dest = 'out_fld',
                        metavar = '', default = wd,
                        help = 'Output folder (.)')

    parser.add_argument('-np', type = int, dest = 'num_p',
                        metavar = 'num_cores', default = 1,
                        help = 'Number of processors to use (1)')

    args = parser.parse_args()

    get_differential_peaks(
        peakfile1 = args.peaks_1,
        peakfile2 = args.peaks_2,
        covfile1 = args.cov_1,
        covfile2 = args.cov_2,
        prefix1 = args.prefix_1,
        prefix2 = args.prefix_2,
        genome_file = args.genome_file,
        window_size = args.window_size,
        stepsize = args.step_size,
        minprobdif = args.min_prob_dif,
        mergedist = args.merge_dist,
        minlen = args.min_len,
        outfld = args.out_fld,
        num_p = args.num_p
    )

if __name__ == "__main__":
    run()
    print(wd)

