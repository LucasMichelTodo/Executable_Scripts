import argparse
import pybedtools as pb
import subprocess as sp
import pathlib as pth
from operator import itemgetter
from collections import defaultdict
from Bio import SeqIO
import os


#### FUNCTIONS ####

def get_extension(filename):
    '''
    Get extension from a filename. It can have up to 2 extensions,
    a filetype extension and optionally a compressing extension (e.g., '.gz')
    '''
    p = pth.Path(filename)
    s = p.suffixes[-2:]
    return(''.join(s))


def sort_gff(gff_file):
    #Transform types to numerical categories to make use of sort.
    types_ord = {
        "protein_coding_gene":1,
        'pseudogene':1,
        "mRNA":2,
        "exon":8,
        "CDS":9,
        "ncRNA_gene":3,
        "rRNA":4,
        "snoRNA":5,
        "snRNA":6,
        "three_prime_UTR":10,
        "tRNA":7,
        "lncRNA":11,
        "asRNA":12
    }

    chr_ord = {
        "Pf3D7_01_v3":1,
        "Pf3D7_02_v3":2,
        "Pf3D7_03_v3":3,
        "Pf3D7_04_v3":4,
        "Pf3D7_05_v3":5,
        "Pf3D7_06_v3":6,
        "Pf3D7_07_v3":7,
        "Pf3D7_08_v3":8,
        "Pf3D7_09_v3":9,
        "Pf3D7_10_v3":10,
        "Pf3D7_11_v3":11,
        "Pf3D7_12_v3":12,
        "Pf3D7_13_v3":13,
        "Pf3D7_14_v3":14,
        "Pf3D7_API_v3":15,
        "Pf3D7_MIT_v3":16
    }

    seqs = []
    with open(gff_file, "r+") as filein:
        for line in filein:
            if line.startswith("#"):
                print(line.strip())
            else:
                seqs.append((line.strip().split("\t")))

    for i in seqs:
        i[0] = chr_ord[i[0]]
        i[2] = types_ord[i[2]]
        i[3] = int(i[3])

    sorted_seqs = sorted(seqs, key=itemgetter(0,3,2))
    for i in sorted_seqs:
        i[0] = [key for key, value in chr_ord.items() if value == i[0]][0]
        i[2] = [key for key, value in types_ord.items() if value == i[2]][0]
        i[3] = str(i[3])

    extension = get_extension(gff_file)
    outgff = gff_file.replace(extension, '_sorted.gff')
    with open(outgff, 'w+') as outfile:
        for i in sorted_seqs:
            outfile.write("\t".join(i)+'\n')

def generate_annotation(ref_gff, types, outfld):

    # Import GFF
    gff = pb.BedTool(ref_gff)

    ## Filter features of interest only
    gene_gff = gff.filter(lambda x: x.fields[2] in types)

    ## Discard Apicoplast
    gene_gff = gene_gff.filter(lambda x: x.chrom != 'Pf3D7_API_v3')

    ## Write output
    p = pth.Path(ref_gff)
    ext = get_extension(ref_gff)
    outname = p.name.replace(ext, '_filtered_only_genes.gff')
    outfld = pth.Path(outfld)
    outfile = pth.Path(outname)
    output = outfld.joinpath(outfile)
    gene_gff.saveas(output)

    ## Sort Gene-GFF
    sort_gff(str(output))

    print(f'Generated: {str(output)}')
    print('And its sorted version.')

def genome_dict(ref_fa):
    # Create genome dict
    genome={}
    records = list(SeqIO.parse(ref_fa, "fasta"))
    genome_file = ref_fa
    for record in records:
        genome[record.id] = (0, len(record.seq))
    return(genome)

def cut_TSS_and_CDS(gff, ref_fa, tss=1000, cds=500):

    genome = genome_dict(ref_fa)

    ref = pb.BedTool(gff)
    current_chrom = ''
    ngenes = len(ref)
    str_bed = ''
    first_in_chrom = False
    last_in_chrom = False

    extension = get_extension(gff)
    suffix = f'_{str(tss)}5p_{str(cds)}CDS.bed'
    outname = gff.replace(extension, suffix)
    with open(outname, 'w+') as outfile:

        for idx, gene in enumerate(ref):

            ## Check Orientation:
            strain = gene.fields[6]

            ## Check if first/last in chromosome
            chrom = gene.chrom

            if current_chrom != chrom:
                first_in_chrom = True
                # print('First in chrom!')
                # print(gene)
                # print('---------------')

            if idx == ngenes-1:
                ## First check if we are in the last gene!
                last_in_chrom = True
            else:
                if ref[idx+1].chrom != chrom:
                    last_in_chrom = True

            ## Set new start and new stop depending on strain:
            if strain == '+':
                newstart = gene.start-tss
                newstop = gene.start+cds
            else:
                newstart = gene.stop-cds
                newstop = gene.stop+tss

            # ## Check overlapp previous gene if +strain or next gene if -strain
            # if strain == '+':
            #     if first_in_chrom:
            #         pass
            #     else:
            #         if newstart < ref[idx-1].stop:
            #             newstart = ref[idx-1].stop+1
            #             # print('Previous gene:')
            #             # print(ref[idx-1])
            #             # print('Current gene:')
            #             # print(gene)
            #             # print('New start-stop:')
            #             # print(newstart, newstop)
            #             # print('--------------------------')
            # else:
            #     if last_in_chrom:
            #         pass
            #     else:
            #         if newstop > ref[idx+1].start:
            #             newstop = ref[idx+1].start-1
            #             # print('Next gene:')
            #             # print(ref[idx+1])
            #             # print('Current gene:')
            #             # print(gene)
            #             # print('New start-stop:')
            #             # print(newstart, newstop)
            #             # print('--------------------------')

            ## Check we dont go < 0
            if newstart < 0: newstart = 0

            ## Check we don't go > chrom length
            if newstop > genome[chrom][1]: newstop = genome[chrom][1]

            ## Check start always < stop
            if newstart >= newstop: newstop = newstart+1

            first_in_chrom = False
            last_in_chrom = False
            current_chrom = chrom

            newline = [gene.chrom, newstart, newstop, gene.fields[8].replace('\t', ' ')]
            newline = [str(x) for x in newline]
            outfile.write('\t'.join(newline)+'\n')
    print(f'Generated: {outname}')

def cross_coverage_with_bed(cov_file, bed_file, outdir):
    '''
    Get coverage from cov_file in regions from bed_file.
    '''
    coverage = pb.BedTool(cov_file)
    bed = pb.BedTool(bed_file)
    name = pth.Path(cov_file).name.replace('.bdg', '')
    bedname = pth.Path(bed_file).name.replace('.bed', '')
    print(f'Calculating {name} coverage over regions in {bedname}.')
    cov = bed.sort().map(coverage, c = 4, o='mean')

    outname = '_crossed_'.join([name, bedname])+'.tsv'
    outfile = pth.Path(outdir).joinpath(pth.Path(outname))
    cov.saveas(outfile)
    print(f'Generated: {str(outfile)}')

def parse_genewise_cov(in_file):
    '''
    Take a '.tsv' genewise coverage file, parse it and return
    a list of lists, one list per line, with [chrom, start, stop,
    coverage, and a dictionary of anotation fileds].
    '''
    parsed_file = []
    with open(in_file, 'r+') as infile:
        for line in infile:
            ll = line.strip().split('\t')
            chrom, start, stop = ll[0:3]
            info = ll[3]
            cov = ll[4]
            fields = info.split(';')
            vals = [f.split('=') for f in fields]
            annot_d = {k:v for k, v in vals}
            parsed_file.append([chrom, start, stop, cov, annot_d])
    return(parsed_file)

#outfile = '/home/lucas/ISGlobal/Projects/Alba/ChIP_Seqs_01_23_edited_genomes/Genewise_Pipe/chip_00160_APC13_genomeDD_q5_sort_noDup_rpkm_normInput_bs10_smth200_pseudo1_crossed_PlasmoDB-61_Pfalciparum3D7_edited_DD_final_filtered_only_genes_10005p_500CDS.tsv'
#x = parse_genewise_cov(outfile)

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

        list_keys = ['GeneTypes', 'Cov_files', 'Sample_names']

        for k, v in pipe_dict.items():
            if k in list_keys:
                outdict[k] = v
            else:
                outdict[k] = v[0]

        return(outdict)
 

def genomewise_coverage(input_file):

    input_d = input_parser(input_file)
    #print(input_d)

    ## Make folders for output
    os.makedirs(input_d['Out_folder'], exist_ok=True)
    outpath = pth.Path(input_d['Out_folder'])
    data_fld = outpath.joinpath('Data')
    os.makedirs(data_fld, exist_ok=True)

    ## Set steps to perform
    make_ref = True if input_d['Run_BuildReference'] == 'yes' else False
    make_bed = True if input_d['Run_BuildBedOverlap'] == 'yes' else False
    cross_cov = True if input_d['Run_CrossCoverage'] == 'yes' else False
    out_table = True if input_d['Run_MakeTable'] == 'yes' else False


    steps = [k for k, v in input_d.items() if k.startswith('Run') and v == 'yes']
    print('Running steps:\n')
    for step in steps:
        print(step)

    if make_ref:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                Generating Gene-wise Reference\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        generate_annotation(
            input_d['Ref_Gff'],
            input_d['GeneTypes'],
            str(data_fld),
        )

    if make_bed:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '             Generate Bed with regions of interest\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        ## Get resulting filename from 'generate annotation'
        p = pth.Path(input_d['Ref_Gff'])
        ext = get_extension(p.name)
        outname = p.name.replace(ext, '_filtered_only_genes.gff')
        outfld = data_fld
        outfile = pth.Path(outname)
        annot_file = str(outfld.joinpath(outfile))

        cut_TSS_and_CDS(
            annot_file,
            input_d['Ref_Fa'],
            tss = int(input_d['Overlap_5p']),
            cds = int(input_d['Overlap_cds'])
        )

    if cross_cov:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                  Cross coverage with bed\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))

        ## Get resulting filename from 'generate annotation'
        p = pth.Path(input_d['Ref_Gff'])
        ext = get_extension(p.name)
        outname = p.name.replace(ext, '_filtered_only_genes.gff')
        outfld = data_fld
        outfile = pth.Path(outname)
        annot_file = str(outfld.joinpath(outfile))
        extension = get_extension(annot_file)
        s5p, scds = input_d['Overlap_5p'], input_d['Overlap_cds']
        suffix = f'_{s5p}5p_{scds}CDS.bed'
        bed_file = annot_file.replace(extension, suffix)

        inpath = pth.Path(input_d['Cov_dir'])
        cov_files = [inpath.joinpath(pth.Path(f)) for f in input_d['Cov_files']]
        for cov_file in cov_files:
            cross_coverage_with_bed(
                cov_file,
                bed_file,
                str(data_fld)
            )

    if out_table:
        print((
            '\n'
            '-------------------------------------------------------------------\n'
            '                    Create output table\n'
            '-------------------------------------------------------------------\n'
            '\n'
        ))
        ## Parse individual files into common table
        cnames = [pth.Path(f).name.replace('.bdg', '') for f in cov_files]
        bname = pth.Path(bed_file).name.replace('.bed', '' )
        outnames = ['_crossed_'.join([cn, bname])+'.tsv' for cn in cnames]
        outfiles = [outfld.joinpath(pth.Path(on)) for on in outnames]
        parsed_outs = [parse_genewise_cov(f) for f in outfiles]

        # parsed outs is a list of lists, each lists represents a line
        # with the order [chrom, start, stop, cov, anot_dict]
        chroms = [x[0] for x in parsed_outs[0]]
        starts = [x[1] for x in parsed_outs[0]]
        stops = [x[2] for x in parsed_outs[0]]
        gids = [x[4]['ID'] if 'ID' in x[4].keys() else 'NA' for x in parsed_outs[0]]
        names = [x[4]['Name'] if 'Name' in x[4].keys() else 'NA' for x in parsed_outs[0]]
        anots = [x[4]['description'] if 'description' in x[4].keys() else 'NA' for x in parsed_outs[0]]
        types = [x[4]['ebi_biotype'] if 'ebi_biotype' in x[4].keys() else 'NA' for x in parsed_outs[0]]
        covs_x = [[x[3] for x in f] for f in parsed_outs]
        covs = [list(x) for x in zip(*covs_x)]

        o_table = 'genomewise_coverage.tsv'
        o_path = pth.Path(input_d['Out_folder'])
        out_file = o_path.joinpath(o_table)
        with open(out_file, 'w+') as outfile:
            cov_cols = ['Coverage_'+x for x in input_d['Sample_names']]
            other_cols = ['Name', 'Annot', 'Ebi_byotype',
                          'Chrom', 'Start', 'Stop']
            cols = ['Gene_id']+cov_cols+other_cols
            outfile.write('\t'.join(cols)+'\n')
            for idx in range(0,len(gids)):
                outline = [gids[idx], *covs[idx], names[idx],
                           anots[idx], types[idx],
                           chroms[idx], starts[idx], stops[idx]]
                outfile.write('\t'.join(outline)+'\n')
        print(f'Final results in: {str(out_file)}')


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

    genomewise_coverage(
        input_file = args.params
    )

if __name__ == "__main__":
    run()
