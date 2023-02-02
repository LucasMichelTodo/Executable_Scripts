import subprocess as sp
import argparse
from statistics import mean


## Create join Bedgraph
def join_bdgs(bdg_list):
    cmd = ['bedtools', 'unionbedg', '-i'] + bdg_list
    union = sp.Popen(cmd, stdout = sp.PIPE)
    str_out = union.stdout.read().decode('utf-8')
    lines = str_out.split('\n')
    return(lines[0:-1])

## Calculate mean coverage
def mean_coverage(bed_list, outname):
    with open(outname, 'w+') as outfile:
        for line in bed_list:
            ll = line.strip().split()
            covs = [float(x) for x in ll[3:]]
            m = mean(covs)
            outfile.write('\t'.join(ll[0:3]+[str(m)])+'\n')

## Main function
def main(bdg_list, outname):
    union = join_bdgs(bdg_list)
    mean_coverage(union, outname)

## Executable Part

def run():
    '''
    Parse command line args and run 'main(args)'.
    '''
    program_description = (
        'Merge Replicates: '
        'Takes a list of bed/bdg files as input and '
        'generates a new bdg file with the mean coverage.'
    )
    parser = argparse.ArgumentParser(description=program_description)

    # Required Arguments
    hline = ('List of space separated \'.bdg\' files.')
    parser.add_argument('-i', nargs='+', type = str,
                        dest = 'bdg_list', required=True,
                        help=hline)

    hline = ('Name for output')
    parser.add_argument('-o', type = str, dest = 'out_name',
                        metavar = 'output_name', required = True,
                        help = hline)

    ## Parse and call main()
    args = parser.parse_args()

    main(
        bdg_list = args.bdg_list,
        outname = args.out_name
    )

if __name__ == "__main__":
    run()

