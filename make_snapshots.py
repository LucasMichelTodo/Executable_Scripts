##############################################################################
##                                                                          ##
##  The idea of this project is to create a python script that takes a      ##
##  list of Gene_ids or a BED file and generates a IGV batch script that    ##
##  takes pictures of each gene/regino automatically.                       ##
##                                                                          ##
##############################################################################

import argparse
import os

wd = os.path.dirname(os.path.realpath(__file__))

def parse_gff_entry_atributes(atributes):
    '''
    Takes the atributes field of a gff entry and parses it.
    Returns a dict with each atribute as key and its value as val.
    '''
    atr_list = atributes.split(';')
    atr_dict = {x.split('=')[0]:x.split('=')[1] for x in atr_list}
    return(atr_dict)

def parse_gff_entry(gff_entry):
    '''
    Takes as input a list with each element of a gff entry (whatever is sepparated by tabs).
    Returns a dict with key-vals for each element.
    '''
    gff_dict = {
        'chrom':gff_entry[0],
        'DB':gff_entry[1],
        'type':gff_entry[2],
        'start':gff_entry[3],
        'stop':gff_entry[4],
        'score':gff_entry[5],
        'strand':gff_entry[6],
        'frame':gff_entry[7],
        'atributes':parse_gff_entry_atributes(gff_entry[8]),
        }
    return(gff_dict)

def parse_gff(gff_file):
    '''
    Reads in a gff file and parses it. It returns a dictionary
    with gene ids as keys and a tupple with chromosome, start and stop as values  
    '''
    gff = {}
    gene_types = ['ncRNA_gene', 'protein_coding_gene', 'pseudogene']
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if line.startswith('#'):
                pass
            else:
                linelist = line.strip().split('\t')
                entry = parse_gff_entry(linelist)
                if entry['type'] in gene_types:
                    gid = entry['atributes']['ID']
                    chrom, start, stop = entry['chrom'], entry['start'], entry['stop']
                    gff[gid] = (chrom, start, stop)
    return(gff)

def parse_gene_ids(ids_file):
    '''
    Reads in a txt file with one gene ID per line and puts the IDs in a list.
    '''
    with open(ids_file, 'r+') as infile:
        gene_ids = [l.strip() for l in infile]

    return(gene_ids)


def main(ref_gff, ids_file, bed_file, files, out_fld, offset, track_h):
    '''
    Main function of the script.
    '''
    geneids = False
    bedfile = False
    os.makedirs(out_fld, exist_ok=True)

    if ref_gff and ids_file:
        gff = parse_gff(ref_gff)
        gene_ids = parse_gene_ids(ids_file)
        geneids = True

    if bed_file: bedfile = True

    ## First common lines
    skeleton = (
        'new\n'
        'expand Gene\n'
        'snapshotDirectory {}\n' ## Posar folder
        'preference IGV.chart.track.height {}\n'
        'setDataRange auto\n'
        .format(out_fld, track_h)
    )

    ## Load files
    loads = ''
    for f in files:
        loads += ('load {}\n' .format(f))

    with open('igv_make_snapshots.bat', 'w+') as outfile:

        outfile.write(skeleton)
        outfile.write(loads)

        ## Snapshot line for each gene id
        if geneids:
            for gid in gene_ids:

                chrom = gff[gid][0]
                start = str(int(gff[gid][1]) - offset)
                stop = str(int(gff[gid][2]) + offset)
                loc = 'goto {}:{}-{}\n' .format(chrom, start, stop)
                comand = 'snapshot\n'
                line = loc+comand
                outfile.write(line)

        ## Snapshot line for each bed line
        if bedfile:
            with open(bed_file, 'r+') as infile:
                for line in infile:
                    ll = line.strip().split()
                    chrom = ll[0]
                    start, stop =  int(ll[1]) - offset, int(ll[2]) + offset
                    loc = 'goto {}:{}-{}\n' .format(chrom, start, stop)
                    comand = 'snapshot\n'
                    line = loc+comand
                    outfile.write(line)

    print('Sanpshot directory: {}' .format(out_fld))
    print(files)


def run():
    '''
    Parse command line args and run 'main'.
    '''
    program_description = ('IGV gene IDs snapshooter: '
                           'Takes as input either a  text file with one gene ID '
                           'per line and a reference GFF_file '
                           'or a BED file with the desired regions, '
                           'and returns an IGV batch script (to run with igv) '
                           'that will take a snapshot of each of the genes/regions.\n\n'
                           'Running this script generates an IGV \'.bat\' file '
                           'which can the be run in IGV using '
                           '\'Tools -> Run Batch Script...\'')

    parser = argparse.ArgumentParser(description=program_description)

    # Required Arguments

    helpline = ('Files to load into IGV. This can be either a single '
                'IGV session file (.xml)(recommended) or a '
                'list of files to load (bam, bed, bdg...) '
                'separated by spaces (file1 file2 file3 ...).')
    parser.add_argument('-f', type = os.path.abspath, dest = 'files',
                        metavar = 'Files', nargs='+', required=True,
                        help=helpline)

    helpline = ('A .txt file with one gene ID per line. '
                'If -ids is set, a reference gff must '
                'be supplied via -gff to get their coordinates.')
    parser.add_argument('-ids', type = str, dest = 'ids_file',
                        metavar = 'IDs_file', required = False,
                        default = False,
                        help = helpline)

    helpline = ('A .gff reference file to get gene coordinates from -ids. '
                'Required only if -ids is set.')
    parser.add_argument('-gff', type = str, dest = 'ref_gff',
                        metavar = 'GFF_file', required = False,
                        default = False,
                        help=helpline)

    parser.add_argument('-bed', type = str, dest = 'bed_file',
                        metavar = 'bed_file', required = False,
                        default = False,
                        help = 'A bed file with desired coordinates for snapshots.')

    # Optional Arguments
    parser.add_argument('-o', '--out_folder', type = os.path.abspath, dest = 'out_fld',
                        metavar = '', default = wd,
                        help = 'Output folder where snapshots will be stored.')

    helpline = ('Number of base-pairs before and after '
        'target gene in the snapshot (3000).')

    parser.add_argument('-offset', type = int, dest = 'offset',
                        metavar = 'Offset', default = 3000,
                        help = helpline)

    parser.add_argument('-th', type = int, dest = 'track_h',
                        metavar = 'Track_height', default = 200,
                        help = 'Track height to be displayed (200).')

    args = parser.parse_args()

    # make one type of input required
    if not args.ids_file and not args.bed_file:
        parser.error("At least one of -ids and -bed required as input.")

    if (args.ids_file and not args.ref_gff) or (args.ref_gff and not args.ids_file):
        parser.error("If -ids is set, -gff is required and vice-versa.")

    main(
        ref_gff = args.ref_gff,
        ids_file = args.ids_file,
        bed_file = args.bed_file,
        files = args.files,
        out_fld = args.out_fld,
        offset = args.offset,
        track_h = args.track_h,
    )


if __name__ == "__main__":
    run()
    print(wd)

