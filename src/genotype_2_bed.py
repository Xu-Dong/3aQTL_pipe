import sys
# transform the 012 genotype matrix into a bed file
# add three columns (chr start end) as the first three columns of genotype matrix file

input_file = sys.argv[1] # the genotype file in 012 format: Genotype_mat.CEU.txt
output_file = sys.argv[2] # the output genotype file in bed format
header_file = sys.argv[3] # store the genotype header
fh = open(input_file,'r')
fho = open(output_file,'w')
fho_header = open(header_file,'w')
header = fh.readline().strip()
print >>fho_header,header
fho_header.close()
i = 0
for line in fh.readlines():
    i += 1
    line = line.strip()
    snp = line.split("\t")[0]
    chrom,pos,ref,alt,release = snp.split("_")
    pos = int(pos)
    print >>fho,"%s\t%d\t%d\t%s" % (chrom,pos-1,pos,line)
    print i

fh.close()
fho.close()
