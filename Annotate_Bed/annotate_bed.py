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


bed = pb.BedTool(bedfile)
ncols_b = len(bed[0].fields)
extracols = ['BedField_' + str(n) for n in range(1, ncols_b-2)]
colnames_b = ['#Bed_Chrom', 'Bed_Start', 'Bed_Stop'] + extracols

gff = pb.BedTool(ref_gff)
gff_cols = ['Gff_Chrom', 'Gff_Source', 'Gff_Feature',
            'Gff_Start', 'Gff_Stop', 'Gff_Score',
            'Gff_Strand', 'Gff_Frame']

cross = bed.intersect(gff, wao=True)
annots = [x.fields[ncols_b+8] for x in cross]
annots = [parse_gff(x) if x is not '.' else {} for x in annots]

## Get all fields present in annotations
## (not all annotations have the same fields)
anot_cols = set(chain.from_iterable([list(x.keys()) for x in annots]))
colnames_g = gff_cols + ['Annot_' + x for x in anot_cols]


with open(outfile, 'w+') as out_file:
    colnames = colnames_b + colnames_g + ['Overlap']
    out_file.write('\t'.join(colnames)+'\n')
    for x in cross:
        info = x.fields[ncols_b+8]
        if info != '.':
            anot_fields = [parse_gff(info).get(y, '.') for y in anot_cols]
        else:
            anot_fields = ['.']*len(anot_cols)
        line = x.fields[:-2] + anot_fields + [x.fields[-1]]
        out_file.write('\t'.join(line)+'\n')


### Calls

wd = '/mnt/Disc4T/Projects/Executable_Scripts/Annotate_Bed/'
os.chdir(wd)
bedfile = 'example_bed.bed'
bedfile = './difpeaks_NF54_+cho_rep2_over_NF54_-cho_rep2_w500_s100_pd0.2_mg500_ml1000_difpeaks.bed'
ref_gff = './PlasmoDB-61_Pfalciparum3D7.gff'
outfile = 'test_output_2.tsv'

