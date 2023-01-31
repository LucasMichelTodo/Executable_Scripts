import pybedtools as pb
import pandas as pd
import numpy as np
import os
from itertools import chain

project_path = "/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/"
os.chdir(project_path)

vep = pb.BedTool("merged_12B_10G_A7_E5_B11_variants_VEPannotated.txt")
gff = pb.BedTool("PlDB-52_Pfalciparum3D7_vep_changetypes.gff.gz")

# Create dict for annotation (from GFF)
gff_gene = gff.filter(lambda x: x[2] in ["gene", "pseudogene"])

def getAnnot(gffentry):

    info = gffentry.fields[8].split(";")
    dinfo = {x.split('=')[0]:x.split('=')[1] for x in info}
    gid = dinfo['ID']
    anot = dinfo['description']
    return([gid, anot])

annot = {}
for entry in gff_gene:
    ga = getAnnot(entry)
    annot[ga[0]] = ga[1]

def getRatioDepth(GF):
    if len(GF) <2:
        rf = np.nan
        alt = np.nan
        ratio = np.nan
        dp = 0
    else:
        rf = int(GF[1].split(",")[0])
        alt = int(GF[1].split(",")[1])
        dp = rf+alt

        if dp == 0:
            ratio = np.nan
        else:
            ratio = round(rf / dp, 2)

    return(rf, alt, ratio, dp)

# Create parsed output

def parse_variant(variant):

    # Parse vcf info
    ref = variant.fields[3]
    alt = variant.fields[4]
    pos = variant.start
    chrom = variant.chrom

    v10G = variant.fields[9].split(":")
    v12B = variant.fields[10].split(":")
    vA7 = variant.fields[11].split(":")
    vB11 = variant.fields[12].split(":")
    vE5 = variant.fields[13].split(":")

    ref_count1, alt_count1, r1, d1 = getRatioDepth(v10G)
    ref_count2, alt_count2, r2, d2 = getRatioDepth(v12B)
    ref_count3, alt_count3, r3, d3 = getRatioDepth(vA7)
    ref_count4, alt_count4, r4, d4 = getRatioDepth(vB11)
    ref_count5, alt_count5, r5, d5 = getRatioDepth(vE5)

    parsed_vcf = [chrom, pos, ref, alt,
                  ref_count1, alt_count1, r1, d1,
                  ref_count2, alt_count2, r2, d2,
                  ref_count3, alt_count3, r3, d3,
                  ref_count4, alt_count4, r4, d4,
                  ref_count5, alt_count5, r5, d5,
                  ]

    # Parse vep info
    info = {}
    for x in variant.fields[7].split(";"):
        feat = x.split("=")
        if len(feat) == 2:
            info[feat[0]] = feat[1]
        else:
            info[feat[0]] = ""

    vep_out = info["CSQ"].split(",")
    effects = [effect.split("|") for effect in vep_out]

    # Add annotation (from GFF)
    for effect in effects:
        gene = effect[4]
        if gene != "":
            gannot = annot[gene]
        else:
            gannot = ""
        effect.append(gannot)

    parsed_variant = [parsed_vcf + effect for effect in effects]

    return(parsed_variant)

# Create DF
colnames = ["Chrom", "Pos", "Ref", "Alt",
            "RefCount_10G", "AltCount_10G", "RefRatio_10G", "depth_10G",
            "RefCount_12B", "AltCount_12B", "RefRatio_12B", "depth_12B",
            "RefCount_A7", "AltCount_A7", "RefRatio_A7", "depth_A7",
            "RefCount_B11", "AltCount_B11", "RefRatio_B11", "depth_B11",
            "RefCount_E5", "AltCount_E5", "RefRatio_E5", "depth_E5",

            "Allele",
            "Consequence",
            "IMPACT",
            "SYMBOL",
            "Gene",
            "Feature_type",
            "Feature",
            "BIOTYPE",
            "EXON",
            "INTRON",
            "HGVSc",
            "HGVSp",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "Existing_variation",
            "DISTANCE",
            "STRAND",
            "FLAGS",
            "SYMBOL_SOURCE",
            "HGNC_ID",
            "SOURCE",
            "PlDB-52_Pfalciparum3D7_vep.gff.gz",

            "Annot"]


parsed = [parse_variant(var) for var in vep]
flat = list(chain.from_iterable(parsed))
var_df = pd.DataFrame.from_records(flat, columns=colnames)

var_df.to_csv("parsed_variants_new.tsv", sep = '\t')
