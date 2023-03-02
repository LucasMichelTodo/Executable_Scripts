import os
import pybedtools as pb
from itertools import chain

def cross_bed_with_gff(cov_file, bed_file, outdir):
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
annots = [x.fields[ncols_bed+8] for x in cross]
annots = [parse_gff(x) if x is not '.' else {} for x in annots]

## Get all fields present in annotations
## (not all annotations have the same fields)
anot_cols = set(chain.from_iterable([list(x.keys()) for x in annots]))
colnames_g = gff_cols + ['Annot_' + x for x in anot_cols]

with open(outfile, 'w+') as out_file:
    out_file.write('\t'.join(colnames_b+colnames_g)+'\n')
    

for x in cross:
    print(x)
    annot_field = x.fields[ncols_bed+8]
    if annot_field != '.':
        annot = parse_gff(annot_field)
        for k, v in annot.items():
            print(f'{k}: {v}')
    print('-------------------------\n')


### Calls

wd = '/mnt/Disc4T/Projects/Executable_Scripts/Annotate_Bed/'
os.chdir(wd)
bedfile = 'example_bed.bed'
bedfile = './difpeaks_NF54_+cho_rep2_over_NF54_-cho_rep2_w500_s100_pd0.2_mg500_ml1000_difpeaks.bed'
ref_gff = './PlasmoDB-61_Pfalciparum3D7.gff'
outfile = 'test_output.tsv'

