import os
import pybedtools as pb
from collections import defaultdict

## FUNTIONS

def filter_by_deletions(del_file, vars_file):

    dd = pb.BedTool(del_file)
    dd_del = dd.filter(lambda x: 'deletion' in x.fields[3])

    dels = defaultdict(list)
    for feat in dd_del:
        dels[feat.chrom ].append((feat.start, feat.stop))

    firstline = True
    out = vars_file.replace('.tsv', '_deletions_filtered.tsv')
    with open(vars_file, 'r+') as infile:
        with open(out, 'w+') as outfile:
            for line in infile:
                if firstline:
                    firstline = False
                    outfile.write(line)
                else:
                    linelist = line.split('\t')
                    chr = linelist[8]
                    pos = int(linelist[9])
                    del_regions = dels[chr]
                    in_del = [pos > _del[0] and pos < _del[1] for _del in del_regions]
                    if not any (in_del):
                            outfile.write('\t'.join(linelist))


## CALLS

wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

dupl_del_dir = '../Paper/Paper_Analysis/Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/'
v_dir = './Parsed_by_Strain/'

dd_files = [    '1.2B_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered.bed',                '10G_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered.bed', 'A7K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered.bed', 'B11_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered.bed', 'E5K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered.bed'
            ]

v_files = [
    '12B_variants_depth_20_refratio_0.5_impactfilter_FALSE_nearest_only.tsv',
    '10G_variants_depth_20_refratio_0.5_impactfilter_FALSE_nearest_only.tsv',
    'A7_variants_depth_20_refratio_0.5_impactfilter_FALSE_nearest_only.tsv',
    'B11_variants_depth_20_refratio_0.5_impactfilter_FALSE_nearest_only.tsv',
    'E5_variants_depth_20_refratio_0.5_impactfilter_FALSE_nearest_only.tsv'
]

for d, v in zip(dd_files, v_files):
    filter_by_deletions(dupl_del_dir+d, v_dir+v)
