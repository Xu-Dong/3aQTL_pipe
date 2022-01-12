import argparse
import os.path

# -- Functions
def load_sample_list(input_sample_list_file):
    sample_list = []
    for line in open(input_sample_list_file,'r'):
        line = line.strip()
        sample_list.append(line)

    return sample_list

def extract_total_reads(input_flagstat_file):
    num_line = 0
    total_reads = '-1'
    #print input_flagstat_file
    for line in open(input_flagstat_file,'r'):
        num_line += 1
        if num_line == 5:
            total_reads = line.strip().split(' ')[0]
            break
    return total_reads


# -- Main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sample_list',help="the file contains all samples")
    parser.add_argument('--path_wig',help="the location of wig files")
    parser.add_argument('--output',help="the final output file with read depth of each sample")

    args = parser.parse_args()

    selected_samples = load_sample_list(os.path.abspath(args.sample_list))
    path_wig = os.path.abspath(args.path_wig)
    
    if path_wig[-1] != "/":
        path_wig += "/"
    else:
        pass
    fho = open(args.output,'w')
    for sample in selected_samples:
        filename_sample = sample + ".flagstat"
        read_depth = extract_total_reads(path_wig + filename_sample)
        wig_location = path_wig + sample + ".wig"
        print >>fho,"%s\t%s" % (wig_location,read_depth)

    fho.close()
