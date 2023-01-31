import os
import pybedtools as pb
from collections import defaultdict

## FUNCTIONS

## Get all 'Consequences'

wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)
raw_vcf = 'parsed_variants.tsv'

firstline = True
consequences = []
cds_consequences = []
with open(raw_vcf, 'r+') as infile:
    for line in infile:
        ll = line.strip().split('\t')
        if firstline:
            colnames = ll
            firstline = False
        else:
            ll = line.strip().split('\t')[1:] #Skip first col (index)
            ld = {k:v for k, v in zip(colnames, ll)}
            consequences.append(ld['Consequence'])
            if ld.get('DISTANCE') == '':
                cds_consequences.append(ld['Consequence'])
colnames
set(consequences)
set(cds_consequences)
cds_vars = set(cds_consequences)
non_cds_vars = set(consequences) - set(cds_consequences)

## Parse variants

def variants_to_dict(var_file):
    variants = {}
    firstline = True
    i = 1
    with open(var_file, 'r+') as infile:
        for line in infile:
            ll = line.strip().split('\t')
            if firstline:
                colnames = ll
                firstline = False
            else:
                ld = {k:[v] for k, v in zip(colnames, ll)}
                var_id = '_'.join(ld['Chrom']+ld['Pos']+ld['Ref']+ld['Alt'])
                if var_id not in variants.keys():
                    variants[var_id] = ld
                    variants[var_id]['Unique_ID'] = f'Variant_{i}'
                    i += 1
                else:
                    ld = {k:v for k, v in zip(colnames, ll)}
                    for k, v in ld.items():
                        variants[var_id][k].append(v)
    return(variants)

## Filter variants


def filter_variants(var_dict, out_file):
    with open(out_file, 'w+') as outfile:

        ## Header
        outfile.write('Var_ID\t'+'\t'.join(colnames)+'\n')

        for k, v in var_dict.items():
            lmask = [vc in cds_vars for vc in v['Consequence']]
            true_idx = [i for i, x in enumerate(lmask) if x]
            if any(lmask):
                for idx in true_idx:
                    row = []
                    for cn in colnames:
                        row.append(v[cn][idx].replace('\"', ''))
                    outfile.write(v['Unique_ID']+'\t'+'\t'.join(row)+'\n')

            elif v['Consequence'] == ['intergenic_variant']:
                row = [v[cn][0].replace('\"', '') for cn in colnames]
                outfile.write(v['Unique_ID']+'\t'+'\t'.join(row)+'\n')

            elif 'upstream_gene_variant' in v['Consequence'] or 'downstream_gene_variant' in v['Consequence']:
                if 'upstream_gene_variant' in v['Consequence']:
                    up_mask = [vc == 'upstream_gene_variant' for vc in v['Consequence']]
                    true_idx = [i for i, x in enumerate(up_mask) if x]
                    to_sort_up = [int(v['DISTANCE'][idx]) for idx in true_idx]
                    for idx, val in enumerate(v['DISTANCE']):
                        if int(val) == min(to_sort_up):
                            up_idx = idx
                    row = [v[cn][up_idx].replace('\"', '') for cn in colnames]
                    outfile.write(v['Unique_ID']+'\t'+'\t'.join(row)+'\n')

                if 'downstream_gene_variant' in v['Consequence']:
                    down_mask = [vc == 'downstream_gene_variant' for vc in v['Consequence']]
                    true_idx = [i for i, x in enumerate(down_mask) if x]
                    to_sort_down = [int(v['DISTANCE'][idx]) for idx in true_idx]
                    for idx, val in enumerate(v['DISTANCE']):
                        if int(val) == min(to_sort_down):
                            down_idx = idx
                    row = [v[cn][down_idx].replace('\"', '') for cn in colnames]
                    outfile.write(v['Unique_ID']+'\t'+'\t'.join(row)+'\n')

## CALLS

wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

var_fld = './Parsed_by_Strain/'
var_fls = [f for f in os.listdir(var_fld) if f.endswith('FALSE.tsv')]

for f in var_fls:
    vars = variants_to_dict(var_fld+f)
    filter_variants(vars, var_fld+f.replace('.tsv', '_nearest_only.tsv'))

# var_file = 'allstrains_variants_depth_20_refratio_0.5_impactfilter_FALSE.tsv'
# filter_variants(
#     variants_to_dict(var_file),
#     var_file.replace('.tsv', '_nearest_only.tsv')
# )
