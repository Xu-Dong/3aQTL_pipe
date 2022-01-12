'''
@ SNP location can be extract directly from processed gt file (embeded in SNP id)
# run on python2 environment
# usage: python extract_SNP_location.py --genotype_bed /xxx/xxx/GEUVADIS.all_chrs.gt012.bed --output /xxx/xxx/output/snp_location.txt
'''
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--genotype_bed',type=str,help="provide the transformed genotype matrix")
parser.add_argument('--output',type=str, default="snp_location.txt",help="specify SNP location file")

args = parser.parse_args()


fh = open(args.genotype_bed,'r')
snp2loc = {} # record the location of snps
snp_order = [] # record all snps to keep order for output

header = fh.readline()
for line in fh.readlines():
    line = line.strip()
    snp = line.split("\t")[0]
    chrom,pos = snp.split("_")[0:2]

    snp2loc[snp] = (chrom,pos)
fh.close()

fho = open(args.output,'w')
print >>fho,"SNP\tChr\tPos"
for snp in snp_order:
    if snp in snp2loc:
        print >>fho,"%s\t%s\t%s" % (snp,snp2loc[snp][0],snp2loc[snp][1])

fho.close()

print "Done!"
