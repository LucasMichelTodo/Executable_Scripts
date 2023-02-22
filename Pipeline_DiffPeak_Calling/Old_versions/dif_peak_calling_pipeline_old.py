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
            'Peaks_1', 'Peaks_2',
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

    #input_file = '/mnt/Disc4T/Projects/Executable_Scripts/Pipeline_DiffPeak_Calling/difpeaks_pipeline_parameters.txt'
    input_d = input_parser(input_file)

    ## Paths
    out_dir = pth.Path(input_d['Out_folder']).resolve()
    peaksdir = pth.Path(input_d['Peaks_dir']).resolve()
    covdir = pth.Path(input_d['Coverage_dir']).resolve()
    os.makedirs(out_dir, exist_ok=True)

    ## Others
    peaks1, peaks2 = input_d['Peaks_1'], input_d['Peaks_2']
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

    peaks = zip(peaks1, peaks2)
    covs = zip(covs1, covs2)
    names = zip(pfxs1, pfxs2)

    for peaks, covs, names in zip(peaks, covs, names):

        # peakfile1 = peaksdir.joinpath(peaks1[0])
        # peakfile2 = peaksdir.joinpath(peaks2[0])
        # covfile1 = covdir.joinpath(covs1[0])
        # covfile2 = covdir.joinpath(covs2[0])
        #prefix1, prefix2 = pfxs1[0], pfxs2[0]

        peakfile1 = peaksdir.joinpath(peaks[0])
        peakfile2 = peaksdir.joinpath(peaks[1])
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

        ## Generate GaussianMixture distribution with 2 components for coverage in peaks
        ## Peaks coverage
        print('Fitting distribution...')

        pf1 = pb.BedTool(peakfile1)
        pf2 = pb.BedTool(peakfile2)
        cf1 = pb.BedTool(covfile1)
        cf2 = pb.BedTool(covfile2)

        x = pf1.map(cf1, c = 4, o='mean')
        for i in range(0,10):
            print(x[i])

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
        
        plot_gaussian_mixture(cov_peaks1, gm1, f'{outfld}/Data/Plots/{prefix1}{suffix}')
        plot_gaussian_mixture(cov_peaks2, gm2, f'{outfld}/Data/Plots/{prefix2}{suffix}')
        plot_components(cov_peaks1, gm1, f'{outfld}/Data/Plots/{prefix1}_components_cdf.png')
        plot_components(cov_peaks2, gm2, f'{outfld}/Data/Plots/{prefix2}_components_cdf.png')

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
        str_beds = str(peakfile1) + ' ' + str(peakfile2)
        outfile = f'{outfld}/Data/Common_Peaks/{prefix1}_{prefix2}_common_peaks.bed'
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

        ## Identify components
        labels1 = gm1.predict(c1s.reshape(-1, 1))
        labels2 = gm2.predict(c2s.reshape(-1, 1))

        comp1_cov1 = np.mean(c1s[labels1 == 0])
        comp2_cov1 = np.mean(c1s[labels1 == 1])

        comp1_cov2 = np.mean(c2s[labels2 == 0])
        comp2_cov2 = np.mean(c2s[labels2 == 1])

        peakcomp1 = 0 if comp1_cov1 > comp2_cov1 else 1
        peakcomp2 = 0 if comp1_cov2 > comp2_cov2 else 1

        print((
            f'Comparing component {peakcomp1+1} from sample {prefix1} '
            f'to component {peakcomp2+1} from sample {prefix2}.'
        ))

        ## Comparing probability of being a peak
        prob1s_pre = gm1.predict_proba(c1s.reshape(-1, 1))
        prob2s_pre = gm2.predict_proba(c2s.reshape(-1, 1))
        prob1s = prob1s_pre[:,peakcomp1]
        prob2s = prob2s_pre[:,peakcomp2]

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

