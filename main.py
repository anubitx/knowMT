from Numt_Tag import numt_chrM_overlapping_reads
from var_ovlp import remove_mitochondrial_noise
from results import get_result
# from tree import contigs_flanking_numts
# from kmer_class import kmer_approach
from alignment_approach import alignment_approach
from cluster import cluster_contigs
import hardcoded_path
import subprocess
import os 
import argparse
import pandas as pd

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-max-ksize", default=35, type=int)
    # parser.add_argument("-min-ksize", default=17, type=int)
    # args = parser.parse_args()
    # n = numt_chrM_overlapping_reads(hardcoded_path.bed_dir, hardcoded_path.bam_dir)
    # n.step3_for_all_files()
    c = remove_mitochondrial_noise()
    c.step3_for_all_files()
    c.get_insert_stats()
    # ref_mmi_cmd = "minimap2 -I 500M -d /home/simarpreet/Desktop/Remote_pc_files/ref.mmi /home/simarpreet/Desktop/Remote_pc_files/original_files/Reference_Genome_GRCh38.p13/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"
    # a = subprocess.run([ref_mmi_cmd], shell = True)
    # out, err = a.communicate()
    # k = kmer_approach(args.max_ksize, args.min_ksize)
    # k = kmer_approach(21,17)
    # k.generate_kmers()
    # k.assemble_masked_reads()
    aln_app = alignment_approach()
    aln_app.write_fastq_files(35)
    aln_app.assemble_masked_reads()
    # flank = contigs_flanking_numts()
    # flank.get_flanking_contigs()  
    cluster = cluster_contigs()
    cluster.get_clustered_contigs()  
    result = get_result()
    result.num_reads_ovlp_MT_nmts()
    result.get_df_from_pickle()
    result.get_result_df()
    print(result.get_all_numts_for_all_samples())

    ##call results.py later on 
    
# a = numt_chrM_overlapping_reads.from_reads_mapped_to_Mt(hardcoded_path.bed_dir, hardcoded_path.Mt_input_bam)


