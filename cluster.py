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
    def __init__(self) -> None:
        self.out = defaultdict(list)
        self.get_minimap_files()
        self.ref_dict = defaultdict(list)
        self.ref_numt_dict()
        self.mt_reads = sorted(os.listdir(hardcoded_path.Mt_reads_dir))
        self.all_data_frame = defaultdict(list)
        # self.ns_ce_dist = []
        # self.ne_cs_dist = []
        # self.SD_range = self.get_insert_info()
        # print(self.SD_range)

    def get_minimap_files(self):
        sub_dirs = sorted(os.listdir(hardcoded_path.Minimap_output))
        for sub_dir in sub_dirs:
            print(sub_dir)
            files = sorted(os.listdir(os.path.join(hardcoded_path.Minimap_output, sub_dir)))
            for file in files:
                self.out[sub_dir].append(file)

    def ref_numt_dict(self):
        for line in open(hardcoded_path.bed_dir_without_chrM):
            fields = line.rstrip().split("\t")
            self.ref_dict[fields[0]].append((fields[1], fields[2], fields[3]))

    def write_to_list(self, distance, query_name, numt_name, tag):
        name_list = []
        name_list.append(distance)
        name_list.append(query_name)
        name_list.append(numt_name)
        name_list.append(tag)
        return name_list

    def dump_df_to_pickle(self):
        if not os.path.exists(os.path.join(hardcoded_path.orignal, 'flanking_contigs_df')):
                os.mkdir(os.path.join(hardcoded_path.orignal, 'flanking_contigs_df'))
        # print(self.all_data_frame)
        for k, v in self.all_data_frame.items():
        # print(self.all_data_frame)
            # idx = 0
            for values in v:
                # print(values)
            # print(v[0])
            # if not os.path.exists(os.path.join(hardcoded_path.orignal, 'flanking_contigs_df', k)):
                # os.mkdir(os.path.join(hardcoded_path.orignal, 'flanking_contigs_df', k))
            # path = os.path.join(hardcoded_path.flank_df, k, f"{k}.df.pkl")

            # idx = 0
            # for file in v:
            #     # print("[df]/df",file)
                # path = os.path.join(hardcoded_path.orignal, 'flanking_contigs_df', k, f"df.{idx+1}.pkl")
                path = os.path.join(hardcoded_path.orignal, 'flanking_contigs_df', f"{k}.flanking_df.pkl")
            #     # final_df.to_pickle(self.all_data_frame[idx],path)
            #     # final_df.to_pickle(path)
                values.to_pickle(path)
                # idx += 1
        
            # print(idx)

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
        
    def get_insert_info(self, k):

        for idx, name in enumerate(self.mt_reads):
            if name.split(".")[0] == k:
                mt_reads = pysam.AlignmentFile(os.path.join(hardcoded_path.Mt_reads_dir, name), 'rb')
                insert = []
                for i in mt_reads:
                    # print(abs(i.template_length) - 200)
                    insert_size = abs(i.template_length) - 200 
                    insert.append(insert_size)      
                mean_insert_size = statistics.mean(insert)
                SD_insert_size = statistics.stdev(insert)
                pvar_insert_size =  statistics.pvariance(insert, mu= mean_insert_size)
                SD_99_low =  mean_insert_size - 3 * SD_insert_size 
                SD_99_upp = mean_insert_size + 3 * SD_insert_size
                return [SD_99_low,SD_99_upp]
        assert(False, "Different file names") 

    def get_clustered_contigs(self):

        # print(ref_dict)
        nc = []
        all_tuples = []
        counter = 0
        for k,v in self.out.items():
            SD_range = self.get_insert_info(k)
            print(SD_range)
            for file in sorted(v): 
                print(file)
                if file.endswith(".bam"):
                    bamfile = pysam.AlignmentFile(os.path.join(hardcoded_path.Minimap_output, k, file), "rb")
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
                    # pd.set_option('display.max_rows', None) 
                    # pd.set_option('display.max_colwidth', None)
                    unique_numts = flank_df['Numt'].unique()
                    print("Total Number of unique reference numts",len(unique_numts))
                    self.all_data_frame[k].append(flank_df)
                    candidate_novel_df = pd.DataFrame(candidate_novel, columns = ["Qname", "Ref_pos", "Ref_name"])
                    candidate_novel_df = candidate_novel_df.sort_values(by = ["Ref_name"])
                    print(candidate_novel_df)
                # self.all_data_frame[k].append(flank_df)
            self.dump_df_to_pickle()
            # break 

cluster = cluster_contigs()
cluster.get_clustered_contigs()
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
    



           # plt.hist(insert)           
        
    #     df = pd.DataFrame(all_tuples, columns=['Chr','Start', 'End', 'Is_Reverse', 'Query_name'])
    #     df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    #     df = df.sort_values(["Start"], ascending= True).groupby(["Chr"])
    #     # a = df.groupby("Chr")
    #     # for i,j in df:
    #         # print(i)
    #         # print(j)
    #     # print(len(all_tuples))

    # # def get_clustered_numts(self):
    #     all_numt_tups = []
    #     for line in open(hardcoded_path.bed_dir_without_chrM):
    #         fields = line.rstrip().split("\t")
    #         # c = list(zip(fields[0], int(fields[1]), int(fields[2]) , fields[3]))
    #         all_numt_tups.append(tuple(fields))
        
    #     numt_df = pd.DataFrame(all_numt_tups, columns=['N_Chr','N_Start', 'N_End', 'Numt_name'])
    #     # print(numt_df)
    #     numt_df[['N_Start', 'N_End']] = numt_df[['N_Start', 'N_End']].astype(int)
    #     numt_df = numt_df.sort_values(["N_Start"], ascending = True).groupby(["N_Chr"])

    #     df['Chr'] == numt_df['N_Chr']
    #     print(df['Query_name'])
        # print(numt_df)
# groupby(["Chr"])
        # for i,j in numt_df:
        #     print(i)
        #     print(j)
        





# cluster = cluster_contigs()
# cluster.get_clustered_contigs()
# read_pickle_file = pd.read_pickle("/home/simarpreet/Desktop/Remote_pc_files/original_files/flanking_contigs_df/ERR1019039/df_1.pkl")
# print(read_pickle_file)
# cluster.get_insert_info()
# cluster.get_clustered_numts()