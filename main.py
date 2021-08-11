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
import shutil

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_file", default=None, type=str)
    parser.add_argument("-numt_bed", default=None, type=str)
    parser.add_argument("-out_dir", default=None, type=str)
    parser.add_argument("-reference", default=None, type=str)
    parser.add_argument("-threads", default=1, type=int)
    parser.add_argument("-ksize", default=35, type=int)
    parser.add_argument("-run_numts_tag", default=False, type=bool)
    # /home/simarpreet/SPAdes-3.12.0-Linux/bin/spades.py
    parser.add_argument("-path_to_spades", default="/home/simarpreet/SPAdes-3.12.0-Linux/bin/spades.py", type=str)    

    args = parser.parse_args()
    input_bam = args.input_file
    input_bed = args.numt_bed
    out_dir = args.out_dir
    reference = args.reference
    threads = args.threads
    ksize = args.ksize
    path_to_spades = args.path_to_spades

    ref_mmi_path = os.path.join(os.getcwd(), 'ref.mmi')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(ref_mmi_path):
        print("................making ref.mmi................")
        # ref_mmi_cmd = f"minimap2 -I 500M -d {ref_mmi_path} /home/simarpreet/Desktop/Remote_pc_files/original_files/Reference_Genome_GRCh38.p13/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"
        ref_mmi_cmd = f"minimap2 -I 500M -d {ref_mmi_path} {args.reference}"
        a = subprocess.run([ref_mmi_cmd], shell = True)
        out, err = a.communicate()
    # Give the path for input arguments mention in Read.me
    # Class 1
    if args.run_numts_tag:
        n = numt_chrM_overlapping_reads(inputbam=input_bam, numtbed=input_bed, 
                                        outputdir=out_dir, threads=threads, reference=reference)
    
    # Class 2
    c = remove_mitochondrial_noise(input_bam=input_bam, output_dir=out_dir, threads=threads, reference= reference)
    c.step3_for_all_files()
    c.get_insert_stats()
    
    # Class 3
    aln_app = alignment_approach(input_bam=input_bam, output_dir=out_dir,
                                 threads=threads, ksize=ksize)
    aln_app.write_fastq_files()
    if args.path_to_spades != None:
        aln_app.assemble_masked_reads(ref_mmi_path, path_to_spades)
    else:
        print("Error::Please provide path to spades.py to run assemble masked reads")

    # Class 4
    cluster = cluster_contigs(input_bam=input_bam, output_dir=out_dir, threads=threads,
                              ksize=ksize, input_bed=input_bed)
    cluster.get_clustered_contigs()  

    # result = get_result(input_bam=input_bam, output_dir=out_dir,
                        # threads=threads, input_bed=input_bed)
    # result.num_reads_ovlp_MT_nmts()
    # result.get_all_numts_for_all_samples()
    # print(result.get_all_numts_for_all_samples())

    ##call results.py later on 

 
# a = numt_chrM_overlapping_reads.from_reads_mapped_to_Mt(hardcoded_path.bed_dir, hardcoded_path.Mt_input_bam)


