import os
import pybedtools as pb
from itertools import chain

def parse_gff(gff_annot):
    '''
    Parse the annotation of a gff file entry. Return it as a dict.
    '''
    fields = gff_annot.split(';')
    vals = [f.split('=') for f in fields]
    annot_d = {k:v for k, v in vals}
    return(annot_d)

def annotate_bed_intersect(bed_file, gff_file, out_file, method = 'closest'):
    bed = pb.BedTool(bed_file)
    ncols_b = len(bed[0].fields)
    extracols = ['BedField_' + str(n) for n in range(1, ncols_b-2)]
    colnames_b = ['#Bed_Chrom', 'Bed_Start', 'Bed_Stop'] + extracols

    gff = pb.BedTool(gff_file)
    gff_cols = ['Gff_Chrom', 'Gff_Source', 'Gff_Feature',
                'Gff_Start', 'Gff_Stop', 'Gff_Score',
                'Gff_Strand', 'Gff_Frame']

    if method == 'intersect':
        cross = bed.intersect(gff, wao=True)
        lastcol = ['Overlap']
    elif method == 'closest':
        cross = bed.closest(gff.sort(), D='b', k=5).saveas()
        lastcol = ['Distance']

    annots = [x.fields[ncols_b+8] for x in cross]
    annots = [parse_gff(x) if x is not '.' else {} for x in annots]

    ## Get all fields present in annotations
    ## (not all annotations have the same fields)
    anot_cols = set(chain.from_iterable([list(x.keys()) for x in annots]))
    colnames_g = gff_cols + ['Annot_' + x for x in anot_cols]


    with open(out_file, 'w+') as outfile:
        colnames = colnames_b + colnames_g + lastcol
        outfile.write('\t'.join(colnames)+'\n')
        for x in cross:
            info = x.fields[ncols_b+8]
            if info != '.':
                anot_fields = [parse_gff(info).get(y, '.') for y in anot_cols]
            else:
                anot_fields = ['.']*len(anot_cols)
            line = x.fields[:-2] + anot_fields + [x.fields[-1]]
            outfile.write('\t'.join(line)+'\n')


### Calls

wd = '/mnt/Disc4T/Projects/Executable_Scripts/Annotate_Bed/'
os.chdir(wd)
bed_file = './difpeaks_NF54_+cho_rep2_over_NF54_-cho_rep2_w500_s100_pd0.2_mg500_ml1000_difpeaks.bed'
#bed_file = 'example_bed.bed'
gff_file = './PlasmoDB-61_Pfalciparum3D7.gff'
out_file = 'test_closest.tsv'

annotate_bed_intersect(bed_file, gff_file, out_file, method='closest')
