from re import template
from typing import DefaultDict
import matplotlib.pyplot as plt
import pysam
import numpy as np
import pandas as pd
from collections import defaultdict
import os
import hardcoded_path
import statistics

class cluster_contigs():
    def __init__(self, input_bam:str, output_dir:str, threads:int, ksize:int, input_bed:str) -> None:
        self.input_bam = input_bam
        self.output_dir = output_dir
        self.threads = threads
        self.input_bed = input_bed
        bambase = os.path.basename(self.input_bam)
        self.ksize = ksize
        self.bambase = bambase
        self.minimap_out = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + f".spadesk{ksize}.minimap.sorted.bam")
        self.mt_reads = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + ".MT.origin.bam")
        self.df_pickle = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + ".flanking_df.pkl")
        self.ref_dict = defaultdict(list)
        self.ref_numt_dict()
        # self.ns_ce_dist = []
        # self.ne_cs_dist = []
        # self.SD_range = self.get_insert_info()
        # print(self.SD_range)

    def ref_numt_dict(self):
        for line in open(self.input_bed):
            fields = line.rstrip().split("\t")
            self.ref_dict[fields[0]].append((fields[1], fields[2], fields[3]))

    def write_to_list(self, distance, query_name, numt_name, tag):
        name_list = []
        name_list.append(distance)
        name_list.append(query_name)
        name_list.append(numt_name)
        name_list.append(tag)
        return name_list

    def get_contig_numt_distance(self, contig_ref_name, contig_start, contig_end, contig_query_name):
        # dist_dict = dict()
        closest_5p = None
        closest_3p = None
        flank_5p= []
        flank_3p = []
        closest_5p_overlap = None
        closest_3p_overlap = None
        overlap_5p = []
        overlap_3p = []
        numt_contains_contig = []
        contig_contains_numt = []
        for key,value in self.ref_dict.items():
            
            if contig_ref_name == key:
                for j in value: 
                    # print(key)
                    ns_cs = int(j[0]) - contig_start
                    ns_ce = int(j[0]) - contig_end
                    ne_cs = int(j[1]) - contig_start
                    ne_ce = int(j[1]) - contig_end
                    # dist_dict[j[2]] = (ns_cs,ns_ce,ne_cs, ne_ce)

                    if ns_cs > 0 and ns_ce > 0 and ne_cs > 0 and ne_ce > 0:
                        if min(ns_cs,ns_ce,ne_cs,ne_ce) == ns_ce:
                            current_closest_5p = ns_ce
                            if closest_5p is not None and current_closest_5p < closest_5p:
                                closest_5p = current_closest_5p
                                flank_5p = self.write_to_list(closest_5p, contig_query_name, j[2], "5p_flank" )
                                # print(flank)
   
                            elif closest_5p is None:
                                closest_5p = current_closest_5p
                                flank_5p = self.write_to_list(closest_5p, contig_query_name, j[2], "5p_flank" )
                                # print(flank)
                
                    elif ns_cs < 0 and ns_ce < 0 and ne_cs < 0 and ne_ce < 0:
                        if max(ns_cs,ns_ce,ne_cs,ne_ce) == ne_cs:
                            current_closest_3p = ne_cs
                            if closest_3p is not None and current_closest_3p > closest_3p:
                                closest_3p = current_closest_3p
                                # print(j[0], j[1])
                                flank_3p = self.write_to_list(abs(closest_3p), contig_query_name, j[2], "3p_flank" )
                                # print(flank_1)
                                
                            elif closest_3p is None:
                                # flank_3p = []
                                closest_3p = current_closest_3p
                                flank_3p = self.write_to_list(abs(closest_3p), contig_query_name, j[2], "3p_flank" )
                                # print(flank_1)
      
                    elif ns_cs > 0 and ns_ce < 0 and ne_cs > 0 and ne_ce > 0:
                        current_closest_5p_overlap = ns_ce
                        if closest_5p_overlap is not None and current_closest_5p_overlap < closest_5p_overlap:
                            closest_5p_overlap = current_closest_5p_overlap
                            overlap_5p = self.write_to_list(closest_5p_overlap, contig_query_name, j[2], "5p_overlap" )
    
                        elif closest_5p_overlap is None:
                            closest_5p_overlap = current_closest_5p_overlap
                            overlap_5p = self.write_to_list(closest_5p_overlap, contig_query_name, j[2], "5p_overlap" )
                
                    elif ns_cs < 0 and ns_ce < 0 and ne_cs > 0 and ne_ce < 0:
                        current_closest_3p_overlap = ne_cs
                        if closest_3p_overlap is not None and current_closest_3p_overlap > closest_3p_overlap:
                            closest_3p_overlap = current_closest_3p_overlap
                            overlap_3p = self.write_to_list( abs(closest_3p_overlap), contig_query_name, j[2], "3p_overlap" )
                    
                        elif closest_3p_overlap is None:
                            closest_3p_overlap = current_closest_3p_overlap
                            overlap_3p = self.write_to_list(abs(closest_3p_overlap), contig_query_name, j[2], "3p_overlap" )

                    #Following condition is based on the assumption that a single contig cannot be contained in 2 numts/ 
                    #unless some reported ref numt events overlap. 
                    elif ns_cs < 0 and ns_ce < 0 and ne_cs > 0 and ne_ce > 0:
                        numt_contains_contig = self.write_to_list(ns_cs, contig_query_name, j[2], "numt_contains_contig" )

                    elif ns_cs > 0 and ns_ce < 0 and ne_cs > 0 and ne_ce < 0:
                        contig_contains_numt = self.write_to_list(ns_cs, contig_query_name, j[2], "contig_contains_numt" )

        return flank_5p, flank_3p, overlap_5p, overlap_3p, numt_contains_contig, contig_contains_numt
        
    def get_insert_info(self):
        mt_reads = pysam.AlignmentFile(self.mt_reads, 'rb', threads=self.threads)
        insert = []
        for i in mt_reads:
            insert_size = abs(i.template_length) - 200 
            insert.append(insert_size)      
        mean_insert_size = statistics.mean(insert)
        SD_insert_size = statistics.stdev(insert)
        SD_99_low =  mean_insert_size - 3 * SD_insert_size 
        SD_99_upp = mean_insert_size + 3 * SD_insert_size
        return [SD_99_low,SD_99_upp]

    def get_clustered_contigs(self):
        nc = []
        all_tuples = []
        counter = 0
        SD_range = self.get_insert_info()
        print(SD_range)
        bamfile = pysam.AlignmentFile(self.minimap_out, 'rb', threads=self.threads)
        # print(os.path.join(hardcoded_path.Minimap_output, k, file))
        distance_5p = []
        distance_3p = []
        num_contigs = []
        final_numt_info = []
        candidate_novel = []
        for contig in bamfile:
            # print(contig)  
            # self.ns_ce_dist = []
            # self.ne_cs_dist = []
            if (not(contig.is_unmapped)):
                num_contigs.append(contig.query_name)
                # print(contig.reference_name, contig.reference_start, contig.reference_end, contig.query_name)
                a = self.get_contig_numt_distance(contig.reference_name, contig.reference_start, contig.reference_end, contig.query_name + "_" + str(counter))
                counter += 1
                if len(a[0]) > 0:
                    # distance_5p.append(a[0][0])
                    if a[0][0] in range(round(SD_range[0]), round(SD_range[1])):
                        final_numt_info.append(a[0])
                if len(a[1]) > 0:
                    # distance_3p.append(a[1][0])
                    # print(a[1])
                    if a[1][0] in range(round(SD_range[0]), round(SD_range[1])):
                        final_numt_info.append(a[1])
                if len(a[2]) > 0:
                    final_numt_info.append(a[2])
                if len(a[3]) > 0:
                    final_numt_info.append(a[3])
                if len(a[4]) > 0:
                    final_numt_info.append(a[4])
                if len(a[5]) > 0:
                    final_numt_info.append(a[5])
                elif len(a[0]) == 0 and len(a[1]) == 0 and len(a[2]) == 0 and len(a[3]) == 0 and len(a[4]) == 0 and len(a[5]) == 0:
                    # if contig.query_sequence != None:
                    # print(contig)
                    candidate_novel.append([contig.query_name, contig.reference_start, contig.reference_name])
        flank_df = pd.DataFrame(final_numt_info, columns=["Distance", "Contig", "Numt", "Overlap/flank"])
        # print(flank_df)
        # pd.set_option('display.max_rows', None) 
        # pd.set_option('display.max_colwidth', None)
        unique_numts = flank_df['Numt'].unique()
        print("Total Number of unique reference numts",len(unique_numts))
        candidate_novel_df = pd.DataFrame(candidate_novel, columns = ["Qname", "Ref_pos", "Ref_name"])
        candidate_novel_df = candidate_novel_df.sort_values(by = ["Ref_name"])
        # print(candidate_novel_df)
        flank_df.to_pickle(self.df_pickle)

# cluster = cluster_contigs()
# cluster.get_clustered_contigs()
            # print(self.all_data_frame)
                # flank_df.to_pickle(os.path.join(os.getcwd(), 'final_df.pkl'))
                # flank3p_df = pd.DataFrame(final_df_3pF, columns=["Distance", "Contig", "Numt", "Overlap/flank"])
                # flank5pO_df = pd.DataFrame(final_df_5pO, columns=["Distance", "Contig", "Numt", "Overlap/flank"])
                # print(flank5p_df)
                # print(flank3p_df)
                # print(flank5pO_df)
                        # print(a[0][0])
                # print(len(num_contigs))
                # print(len(distance_5p))

                # print(min(distance_3p), max(distance_3p))
                # # print(sum(distance_5p) / len(distance_5p))
                # mean_insert_size = statistics.mean(distance_3p)
                # SD_insert_size = statistics.stdev(distance_3p)
                # print(mean_insert_size , "Mean")
                # print("SD", SD_insert_size)
                
                # print(SD_insert_size - mean_insert_size)
                # print(mean_insert_size + SD_insert_size)
                # print(3 * SD_insert_size - mean_insert_size)
                # print(mean_insert_size + 3 * SD_insert_size)
                # print(set(distance_5p))
                # plt.hist(distance_5p)
                # plt.hist(distance_5p, range = (0, 10000))
                # plt.title("Distribution of distance to the numt from the closest 5p contig")
                # plt.xlabel("Distance in base pairs")
                # plt.ylabel("Number of Contigs")
                # plt.savefig("5p_distance_from_numts.png")


                # print(distance_5p)
            
                    # break
    
# cluster = cluster_contigs()
# cluster.get_clustered_contigs()
# read_pickle_file = pd.read_pickle("/home/simarpreet/Desktop/Remote_pc_files/original_files/flanking_contigs_df/ERR1019039/df_1.pkl")
# print(read_pickle_file)
# cluster.get_insert_info()
# cluster.get_clustered_numts()