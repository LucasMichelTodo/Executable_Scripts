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
from collections import defaultdict
import pathlib as pth
from pathlib import Path
from scipy.stats import multivariate_normal


## Functions

def create_window_ref(window_size, genome_file, step):
    bed = pb.BedTool()
    win_ref = bed.window_maker(w=window_size, g = genome_file, s = step)
    return(win_ref)

def plot_gaussian_mixture(cov_vector, gm, plotname):

    # Plot the histogram.
    plt.hist(cov_vector, bins=250, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    logprob = gm.score_samples(x.reshape(-1, 1))
    pdf = np.exp(logprob)
    plt.plot(x, pdf, 'k', linewidth=2, label = 'Join PDF')

    # derive cumulative distribution function (cdf)
    cumsum = np.cumsum(pdf)
    # scale as a probability distribution
    cdf = cumsum / np.max(cumsum)
    # plot cdf, scaled to the y limits of the above plot
    xmin, xmax, ymin, ymax = plt.axis()
    plt.plot(x, cdf * ymax, '-b', label='Join CDF') # multiplying by ymax sets it to plot scale

    ## Get each distribution and it's associated probabilities
    multi_var_distr = [multivariate_normal(mu, sigma) for (mu, sigma) in zip(gm.means_, gm.covariances_)]
    cdf1 = multi_var_distr[0].cdf(x)
    cdf2 = multi_var_distr[1].cdf(x)
    # multiplying by ymax sets it to plot scale
    plt.plot(x, cdf1* ymax, ':r', linewidth=2, label = 'Component 1 CDF')
    # multiplying by ymax sets it to plot scale
    plt.plot(x, cdf2* ymax, ':g', linewidth=2, label = 'Component 2 CDF')

    plt.legend(loc='upper left')
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
    plt.plot(x, cdf1, 'r', linewidth=2, label="Component 1")
    plt.plot(x, cdf2, 'g', linewidth=2, label="Component 2")
    plt.legend(loc='upper left')
    plt.savefig(plotname)
    plt.clf()

def get_peak_component(gm, cov):
    '''
    Take a gm model and a coverage vector and return components
    by coverage mean order (so first returned component
    always corresponds to peaks).
    '''
    labels = gm.predict(cov.reshape(-1,1))
    cmp0cov = np.mean(cov[labels == 0])
    cmp1cov = np.mean(cov[labels == 1])

    peakcmp = 0 if cmp0cov > cmp1cov else 1
    bkgdcmp = 1 - peakcmp
    return([peakcmp, bkgdcmp])

def get_individual_peaks(covfile, gm, mergedist, minlen, outfld, name):

    cov_file = pb.BedTool(covfile)
    cov = np.array([float(f.fields[3]) for f in cov_file])
    probs = gm.predict_proba(cov.reshape(-1,1))
    peakcomp, bkgdcomp = get_peak_component(gm, cov)
    mask = probs[:,bkgdcomp] < probs[:,peakcomp]

    outbed = Path(outfld).joinpath(f'{name}_individual_peaks.bed')
    with open(outbed, 'w+') as outfile:
        for feat, mask in zip(cov_file, mask):
            if mask:
                outfile.write(str(feat))

    outbed_raw = pb.BedTool(outbed)
    out_pre = outbed_raw.sort().merge(d=mergedist)
    out = out_pre.filter(lambda f: f.stop - f.start > minlen).saveas(outbed)

def plot_scores(score_vect, chroms, starts, stops, outfile):
        with open(outfile, 'w+') as scorebed:
            for i in range(0, len(score_vect)):
                ll = [chroms[i], starts[i], stops[i], score_vect[i]]
                scorebed.write('\t'.join([str(x) for x in ll])+'\n')

def add_peak_names_reorder_cols(bed, outname, header):
    names = ['Peak_']*len(bed)
    idxs = range(1,len(bed)+1)
    peaknames = [n+str(i) for n,i in zip(names, idxs)]

    with open(outname, 'w+') as outfile:
        outfile.write(header)
        for feat, name in zip(bed, peaknames):
            l = feat.fields
            outlist = l[0:3] + [name]+ [l[-1]] + l[3:-1]
            outfile.write('\t'.join(outlist)+'\n')

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

        list_keys = [
            'Coverages_1', 'Coverages_2',
            'Names_1', 'Names_2'
        ]

        for k, v in pipe_dict.items():
            if k in list_keys:
                outdict[k] = v
            else:
                outdict[k] = v[0]

        return(outdict)


## Diffpeakcalling algortihm
def get_differential_peaks(input_file):

    input_d = input_parser(input_file)
    print('\nProgram parameters:\n')
    for k, v in input_d.items():
        print(f'{k}:')
        print(f'{v}\n')


    ## Paths
    out_dir = Path(input_d['Out_folder']).resolve()
    covdir = Path(input_d['Coverage_dir']).resolve()
    os.makedirs(out_dir, exist_ok=True)

    ## Others
    covs1, covs2 = input_d['Coverages_1'], input_d['Coverages_2']
    pfxs1, pfxs2 = input_d['Names_1'], input_d['Names_2']
    genome_file = input_d['Genome']
    window_size = int(input_d['Win_Size'])
    stepsize = int(input_d['Step'])
    minprobdif = float(input_d['MinProbDif'])
    mergedist = int(input_d['Merge'])
    minlen = int(input_d['MinLen'])
    num_p = int(input_d['Threads'])

    start = time.time()

    ## Greet User
    print((
        '\n'
        '###################################################################\n'
        '###   Running Custom differential Peak-Calling by Lucas M.T.    ###\n'
        '###################################################################\n'
        '\n'
    ))

    # Run Algorithm by Pairs

    covs = zip(covs1, covs2)
    names = zip(pfxs1, pfxs2)

    for covs, names in zip(covs, names):

        #covfile1 = covdir.joinpath(covs1[0])
        #covfile2 = covdir.joinpath(covs2[0])
        #prefix1, prefix2 = pfxs1[0], pfxs2[0]

        covfile1 = covdir.joinpath(covs[0])
        covfile2 = covdir.joinpath(covs[1])
        prefix1, prefix2 = names[0], names[1]

        print(f'\nComparing samples: {prefix1} vs {prefix2}\n')
        outfld = out_dir.joinpath(prefix1+'_vs_'+prefix2)

        ## Create folders for output
        os.makedirs(outfld, exist_ok = True)
        os.makedirs(outfld.joinpath('Data/'), exist_ok = True)
        os.makedirs(outfld.joinpath('Data/Common_Coverage/'), exist_ok = True)
        os.makedirs(outfld.joinpath('Data/Common_Peaks/'), exist_ok = True)
        os.makedirs(outfld.joinpath('Data/Plots/'), exist_ok = True)
        os.makedirs(outfld.joinpath('Data/Window_Coverage/'), exist_ok = True)
        os.makedirs(outfld.joinpath('Data/PreFilter_Difpeaks/'), exist_ok = True)

        ## Generate GaussianMixture distribution with 2 components for coverage

        print('Fitting distribution...')

        cf1 = pb.BedTool(covfile1)
        cf2 = pb.BedTool(covfile2)
        cov_1 = np.array([float(f.fields[3]) for f in cf1])
        cov_2 = np.array([float(f.fields[3]) for f in cf2])

        ## Fitting distribution: Gaussian Mixture

        ## All sklearn estimators use as input a 2D array,
        ## with samples as rows and features as columns.
        ## Our data-set consists in one sample per bed feature
        ## with only one feature (coverage), so we have to reshape it
        ## into a 2D array with just one column.
        ## -1 in reshape just means whatever it takes to make it work!
        ## So (-1, 1) means: any number of rows and 1 column

        gm1 = GaussianMixture(n_components=2).fit(cov_1.reshape(-1, 1))
        gm2 = GaussianMixture(n_components=2).fit(cov_2.reshape(-1, 1))

        suffix = '_peak_coverage_fit.png'

        # Plotting Fitted Distribution

        plot_gaussian_mixture(cov_1, gm1, f'{outfld}/Data/Plots/{prefix1}{suffix}')
        plot_gaussian_mixture(cov_2, gm2, f'{outfld}/Data/Plots/{prefix2}{suffix}')
        plot_components(cov_1, gm1, f'{outfld}/Data/Plots/{prefix1}_components_cdf.png')
        plot_components(cov_2, gm2, f'{outfld}/Data/Plots/{prefix2}_components_cdf.png')

        ## Identify components
        ## (one component represents peak regions and the other background)
        peakcomp1, bkgdcomp1 = get_peak_component(gm1, cov_1)
        peakcomp2, bkgdcomp2 = get_peak_component(gm2, cov_2)

        ## Generating individual peaks
        get_individual_peaks(covfile1, gm1, 100, 100, outfld, prefix1)
        get_individual_peaks(covfile2, gm2, 100, 100, outfld, prefix2)

        ## Create window reference
        print('Creating windowed coverage...')
        suffix = '_window_coverage.bdg'
        win_ref = create_window_ref(window_size, genome_file, stepsize)
        win_cov1 = win_ref.map(cf1, c=4, o='mean').saveas(f'{outfld}/Data/Window_Coverage/{prefix1}{suffix}')
        win_cov2 = win_ref.map(cf2, c=4, o='mean').saveas(f'{outfld}/Data/Window_Coverage/{prefix2}{suffix}')

        ## Create join coverage bdg
        print('Creating join coverage bed...')
        outfile = outfld.joinpath(f'Data/Common_Coverage/{prefix1}_{prefix2}_common_coverage.bdg')
        cmd = ['bedtools', 'unionbedg', '-i',
               f'{outfld}/Data/Window_Coverage/{prefix1}{suffix}',
               f'{outfld}/Data/Window_Coverage/{prefix2}{suffix}',
               '>', str(outfile)]
        sp.call(' '.join(cmd), shell = True)
        union_bed_pre = pb.BedTool(outfile)

        ## Filter in case there are conflicting regions in genome
        ## (apicoplast with different name i.e.). Conflicting
        ## regions will be ignored though!
        union_bed = union_bed_pre.filter(lambda b: b.fields[3] != '.' and b.fields[4] != '.')

        ## Create common peaks bed
        print('Creating common peaks bed...')

        files_to_cross = []
        peakfile1 = outfld.joinpath(f'{prefix1}_individual_peaks.bed')
        peakfile2 = outfld.joinpath(f'{prefix2}_individual_peaks.bed')
        str_beds = str(peakfile1) + ' ' + str(peakfile2)
        outfile = f'{outfld}/Data/Common_Peaks/{prefix1}_{prefix2}_common_peaks.bed'
        cmd = 'awk \'{print}\' '+f'{str_beds} > {outfile}'
        sp.call(cmd, shell = True)
        common_peaks = pb.BedTool(outfile).sort()

        ## Subset union_bed to common_peaks
        common_cov_common_peaks = union_bed.intersect(common_peaks)

        ## Call differential peaks

        print('Calling differential peaks (this might take a while)...')

        c1s = np.array([float(x.fields[3]) for x in common_cov_common_peaks])
        c2s = np.array([float(x.fields[4]) for x in common_cov_common_peaks])

        print((
            f'Comparing component {peakcomp1+1} from sample {prefix1} '
            f'to component {peakcomp2+1} from sample {prefix2}.'
        ))

        # Calculate CDF, for each interval for each component

        multi_var_distr1 = [multivariate_normal(mu, sigma) for (mu, sigma) in zip(gm1.means_, gm1.covariances_)]
        cdf1_peakcomp = multi_var_distr1[peakcomp1].cdf(c1s)
        cdf1_bkgdcomp = multi_var_distr1[bkgdcomp1].cdf(c1s)

        multi_var_distr2 = [multivariate_normal(mu, sigma) for (mu, sigma) in zip(gm2.means_, gm1.covariances_)]
        cdf2_peakcomp = multi_var_distr2[peakcomp2].cdf(c2s)
        cdf2_bkgdcomp = multi_var_distr2[bkgdcomp2].cdf(c2s)

        cdf_peak_dif = cdf1_peakcomp - cdf2_peakcomp
        cdf_bkgd_dif = cdf1_bkgdcomp - cdf2_bkgdcomp

        # Plot peaks-component cdf difference distribution
        plotname = f'{outfld}/Data/Plots/cdf_peakcomp_dif.png'
        plt.hist(cdf_peak_dif, bins=25, density=True, alpha=0.6, color='r')
        plt.savefig(plotname)
        plt.clf()

        # Plot background-component cdf difference distribution
        plotname = f'{outfld}/Data/Plots/cdf_bkgdcomp_dif.png'
        plt.hist(cdf_bkgd_dif, bins=25, density=True, alpha=0.6, color='r')
        plt.savefig(plotname)
        plt.clf()

        ## Getting Peaks

        peaks1over2_peaks = cdf_peak_dif >= minprobdif
        peaks1over2_bkgd = cdf_bkgd_dif >= minprobdif
        peaks2over1_peaks = cdf_peak_dif <= -minprobdif
        peaks2over1_bkgd = cdf_bkgd_dif <= -minprobdif

        peaks1over2 = [any(b) for b in zip(peaks1over2_peaks, peaks1over2_bkgd)]
        peaks2over1 = [any(b) for b in zip(peaks2over1_peaks, peaks2over1_bkgd)]

        chroms = [x.chrom for x in common_cov_common_peaks]
        starts = [x.start for x in common_cov_common_peaks]
        stops = [x.stop for x in common_cov_common_peaks]

        ## Create scoring function bdg
        outplot = f'{outfld}/Data/Plots/cdf_peakcomp_dif.bdg'
        plot_scores(cdf_peak_dif, chroms, starts, stops, outplot)
        outplot = f'{outfld}/Data/Plots/cdf_bkgdcomp_dif.bdg'
        plot_scores(cdf_bkgd_dif, chroms, starts, stops, outplot)

        ## Generate final list of Differential Peaks
        pd.DataFrame(common_cov_common_peaks)
        df = pd.read_table(common_cov_common_peaks.fn, names=['#chrom', 'start', 'stop', 'cov1', 'cov2'])
        df['Bkgd_CDF_dif'] = cdf_bkgd_dif
        df['Peak_CDF_dif'] = cdf_peak_dif
        ## Get max (in positive or negative)
        row_max = df[['Bkgd_CDF_dif', 'Peak_CDF_dif']].abs().max(axis=1)
        df['Max_CDF_dif'] = df[['Bkgd_CDF_dif', 'Peak_CDF_dif']].max(axis=1).mask(lambda x: x < row_max, -row_max)

        rawout1 = (f'{outfld}/Data/PreFilter_Difpeaks/'
                   f'difpeaks_{prefix1}over{prefix2}'
                   f'_w{window_size}_s{stepsize}_pd{minprobdif}.bed')

        rawout2 = (f'{outfld}/Data/PreFilter_Difpeaks/'
                   f'difpeaks_{prefix2}over{prefix1}'
                   f'_w{window_size}_s{stepsize}_pd{minprobdif}.bed')

        df[peaks1over2].to_csv(rawout1, sep='\t', index=False)
        df[peaks2over1].to_csv(rawout2, sep='\t', index=False)

        rawbed1 = pb.BedTool(rawout1)
        rawbed2 = pb.BedTool(rawout2)

        ## Merge peaks and filter by length
        out1 = (f'{outfld}/difpeaks_{prefix1}_over_{prefix2}'
                f'_w{window_size}_s{stepsize}_pd{minprobdif}'
                f'_mg{mergedist}_ml{minlen}'
                '_difpeaks.bed')

        out2 = (f'{outfld}/difpeaks_{prefix2}_over_{prefix1}'
                f'_w{window_size}_s{stepsize}_pd{minprobdif}'
                f'_mg{mergedist}_ml{minlen}'
                '_difpeaks.bed')

        funcs = ['mean']*5
        out1_pre = rawbed1.sort().merge(d=mergedist, c =[4,5,6,7,8], o=funcs)
        bed1 = pb.BedTool(out1_pre.filter(lambda f: f.stop - f.start > minlen).saveas())
        out2_pre = rawbed2.sort().merge(d=mergedist, c =[4,5,6,7,8], o=funcs)
        bed2 = pb.BedTool(out2_pre.filter(lambda f: f.stop - f.start > minlen).saveas())

        ## Add header, peaknames and reorder columns
        header = (
            f'# Differential peaks for {prefix1} vs {prefix2}\n'
            '# This file was generated to be run in a browser like IGV.\n'
            '# Columns correspond to: chr, start, stop, peak_id, max_cdf_dif, '
            'coverage sample1, coverage sample2, bkgd_comp_cdf_dif, '
            'peak_cmp_cdf_dif, max_cdf_dif.\n'
        )
        add_peak_names_reorder_cols(bed1, out1, header)
        add_peak_names_reorder_cols(bed2, out2, heade
                                    )

    end = time.time()
    print('Finished! Elapsed time:')
    print(end - start)

## Parse Arguments and run program

wd = os.path.dirname(os.path.realpath(__file__))

def run():
    '''
    Parse parameters file and run 'get_differential_peaks(args)'.
    '''
    program_description = ('Differential Peaks Caller:'
                           'For each sample, takes as input a \'.narrowPeak\' file '
                           '(from a MACS2 callpeak call) '
                           'and a bed/bedgraph coverage file '
                           '(generated with DeepTools e.g.) '
                           'and returns a list of differential peaks.')

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

    get_differential_peaks(
        input_file = args.params
    )

if __name__ == "__main__":
    run()
    print(wd)

