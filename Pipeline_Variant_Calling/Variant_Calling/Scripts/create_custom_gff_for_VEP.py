import pybedtools as pb
import subprocess as sp
import os

outfld = "/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/"
os.chdir(outfld)

# Filter out lines that do not correspond to genes
gff = pb.BedTool("../Data/PlasmoDB-52_Pfalciparum3D7.gff")

types = [feat.fields[2] for feat in gff]
set(types)


# Add biotype info to transcripts (it is a VEP requisite)
biotypes = {
    "mRNA": "protein_coding",
    "ncRNA": "ncRNA",
    "rRNA": "rRNA",
    "snoRNA": "snoRNA",
    "snRNA": "snRNA",
    "tRNA": "tRNA",
    "five_prime_UTR": "five_prime_UTR",
    "three_prime_UTR": "three_prime_UTR",
    "pseudogenic_transcript": "pseudogenic_transcript"
}

def add_biotype(feature):
    if feature.fields[2] in biotypes.keys():
        feature.attrs["biotype"] = biotypes[feature.fields[2]]
    return(feature)

added_biotype = gff.each(add_biotype)

# Sort GFF

added_biotype.sort().saveas("PlDB-52_Pfalciparum3D7_vep.gff")

# Change type from "protein_coding_gene" to "gene"

types = ['ncRNA_gene', 'protein_coding_gene']

with open("PlDB-52_Pfalciparum3D7_vep.gff", 'r+') as infile:
    with open("PlDB-52_Pfalciparum3D7_vep_changetypes.gff", 'w+') as outfile:
        for line in infile:
            linelist = line.strip().split('\t')
            if linelist[2] in types:
                linelist[2] = 'gene'
            outfile.write('\t'.join(linelist)+'\n')

# Compress GFF
cmd = "bgzip PlDB-52_Pfalciparum3D7_vep_changetypes.gff"
sp.call(cmd, shell=True)

# Tabix GFF
cmd = "tabix -p gff PlDB-52_Pfalciparum3D7_vep_changetypes.gff.gz"
sp.call(cmd, shell=True)
