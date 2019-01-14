from HTSeq import GenomicInterval, GenomicArray, BED_Reader
from sys import argv,exit
import peakFeatures
usage = """
python calPeakFeatures.py bedfile bdgfile outfile
#1. bed file
#2. bedgraph file
"""
if len(argv) != 4:
    print(usage)
    exit()

bedfile = argv[1]
bdgfile = argv[2]
outfile = argv[3]
#chroms = list(map(lambda x:"chr"+str(x), range(1,23))) + ["chrX","chrY","chrM"]

def build_genome_array(bdgfile):
    #genome_array = GenomicArray(chroms,stranded=False,typecode="d")
    genome_array = GenomicArray("auto", stranded=False, typecode="d")
    with open(bdgfile) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            chrom,start,end,value = line.strip().split()
            iv = GenomicInterval(chrom,int(start),int(end),".")
            genome_array[iv] += float(value)
    return genome_array

genome_array = build_genome_array(bdgfile)
peaks = BED_Reader(bedfile)
with open(outfile,"w") as fout:
    #fout.write("peakName" + "\t" + "\t".join(peakFeatures.peakFeaturesHeader()) + "\n")
    for i in peaks:
        values = genome_array[i.iv].values()
        values = list(map(int,values))

        features = peakFeatures.peakFeatures(values)
        #fout.write(i.name + "\t" + "\t".join(map(str,features)) + "\n")
        fout.write("\t".join(map(str, features)) + "\n")

