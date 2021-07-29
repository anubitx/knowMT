from Numt_Tag import numt_chrM_overlapping_reads
# from var_ovlp import remove_mitochondrial_noise
# from tree import contigs_flanking_numts
# from kmer_class import kmer_approach
import hardcoded_path
import argparse
import shutil, os

def copy_bam_files(bam_folderpath):
    destination = hardcoded_path.bam_dir
    for i in os.listdir(bam_folderpath):
        shutil.copy(i, destination)  
    # print(os.listdir(destination))
        # print(i)

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-max-ksize", default=35, type=int)
    parser.add_argument("-threads", default=1, type=int)
    parser.add_argument("-min-ksize", default=17, type=int)
    parser.add_argument("-numtbed", default=None, type=str)
    parser.add_argument("-inputbam", default=None, type=str)
    parser.add_argument("-outputdir", default=None, type=str)
    parser.add_argument("-reference", default=None, type=str)
    args = parser.parse_args()
    # if args.bam_folderpath is not None:
    #     copy_bam_files(args.bam_folderpath)
    # if args.bed_filename is not None:
    #     shutil.copy(args.bed_filename, hardcoded_path.bed_dir) 

    a = numt_chrM_overlapping_reads(args.numtbed, args.inputbam, args.outputdir, args.threads, args.reference)
    # a = numt_chrM_overlapping_reads(args.bed_filename, args.bam_folderpath)
    # copy_bam_files("original_files/Mt_input_bam")
    # a = numt_chrM_overlapping_reads(args.bed_file_name, args.bam_file_path)
    
    # a.step3_for_all_files()
    # a = numt_chrM_overlapping_reads.from_reads_mapped_to_Mt(hardcoded_path.bed_dir, hardcoded_path.Mt_input_bam)
    # # a.generate_numt_tag()
    # c = remove_mitochondrial_noise(hardcoded_path.step3_output_dir)
    # c.get_insert_stats()
    # k = kmer_approach(args.max_ksize, args.min_ksize)
    # k.generate_kmers()
    # k.assemble_masked_reads()
    # flank = contigs_flanking_numts()
    # flank.get_flanking_contigs()    