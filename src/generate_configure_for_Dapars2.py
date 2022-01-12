import argparse
import os.path

# -- Functions
def extract_all_wigs(input_seq_depth_file):
    wig_files = []
    for line in open(input_seq_depth_file,'r'):
        line = line.strip()
        w = line.split("\t")
        wig_files.append(w[0])

    return wig_files


# -- Main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--annotation_3utr',help="the location of reference 3'UTR bed file")
    parser.add_argument('--wigFile_depth',help="the index file contains all wig files and read depth")
    parser.add_argument('--coverage_threshold',help="specify the threshold of coverage")
    parser.add_argument('--threads',help="specify the number of threads used")
    parser.add_argument('--out_dir',help="specify the directory of dapars2 output,do not use absolute path")
    parser.add_argument('--out_prefix',help="specify the name of result file of dapars2")
    parser.add_argument('--out_config_name',help="specify configure file name")

    args = parser.parse_args()

    configure_file_name = args.out_config_name
    fho = open(configure_file_name,'w')
    # print Annotated_3UTR
    print >>fho,"# Specify the reference of 3'UTR region"
    print >>fho,"\nAnnotated_3UTR=" + os.path.abspath(args.annotation_3utr)

    # print wig files
    all_wig_files = extract_all_wigs(os.path.abspath(args.wigFile_depth))
    print >>fho,"\n# A comma separated list of wig files of all samples"
    print >>fho,"\nAligned_Wig_files=" + ",".join(all_wig_files)
    
    # specify Output_directory and Output_result_file
    print >>fho,"\nOutput_directory=" + args.out_dir
    print >>fho,"\nOutput_result_file=" + args.out_prefix

    # specify Coverage_threshold
    print >>fho,"\n# Specify Coverage threshold"
    print >>fho,"\nCoverage_threshold=" + args.coverage_threshold

    # specify the Num_Threads to process the analysis
    print >>fho,"\n# Specify the numbero of threads to process the analysis"
    print >>fho,"\nNum_Threads=" + args.threads

    # Provide sequencing_depth_file for normalization
    print >>fho,"\n# Provide sequencing depth file for normalization"
    print >>fho,"\nsequencing_depth_file=" + os.path.abspath(args.wigFile_depth)
    fho.close()
    

