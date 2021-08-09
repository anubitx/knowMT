import os

from numpy.lib.function_base import i0
import hardcoded_path
import pysam
import pandas as pd
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import re

class get_result():
    def __init__(self):
        self.outfiles = sorted(os.listdir(hardcoded_path.out_dir))
        self.mt_reads = sorted(os.listdir(hardcoded_path.Mt_reads_dir))
        self.numt_reads = sorted(os.listdir(hardcoded_path.numt_reads))
        self.flank_df = sorted(os.listdir(hardcoded_path.flank_df))
        # self.bed_dir = os.listdir(hardcoded_path.bed_dir_without_chrM)
        self.numts_1read = []
        self.all_files_total_nmts = []
        

    def num_reads_ovlp_MT_nmts(self):
        all_files = []
        for idx, name in enumerate(self.outfiles):
            num_ovlp_reads = [] 
            count_nmts = []
            actual_numts = set()
            samfile = pysam.AlignmentFile(os.path.join(hardcoded_path.out_dir,name),'rb')
            # print(name)
            for aln in samfile:
                num_ovlp_reads.append(aln.query_name)
                tag = aln.get_tag("Ni")
                count_nmts.append(tag.split(",")[2])
                actual_numts.add(tag.split(",")[2])
                # actual_numts.add(tag.split(",")[7])
            print("actual num of numts",len(actual_numts) - 3)
            self.all_files_total_nmts.append(len(actual_numts) - 3)

            # print(Counter(count_nmts))
            read_count = Counter(count_nmts)
            res = [k for k,v in read_count.items() if v == 1]
            print(res)
            self.numts_1read.append(res)
            counts_hist = [v for k,v in read_count.items()]
            plt.hist(counts_hist, range = (0,10))
            plt.title("Distribution of Numts by number of supporting reads")
            plt.xlabel("No. of supporting reads")
            plt.xticks(range(11))
            plt.ylabel("No. of Numts")
            plt.savefig(f"Num_of_supporting_reads{idx}.png")
            plt.close()
            all_files.append((name.split(".")[0], len(num_ovlp_reads)))
            # print( name.split(".")[0],"No. of reads that overlap nmts/chrM" ,len(num_ovlp_reads))
        # print(all_files)
        
        
        # print(all_files_total_nmts)
        plt.bar(range(len(self.all_files_total_nmts)), self.all_files_total_nmts, color = 'g')
        for index,data in enumerate(self.all_files_total_nmts):
                plt.text(x=index , y =data , s=f"{data}" , fontdict=dict(fontsize=10), ha = 'center')
        plt.title("Total number of numts with atleast 1 supporting read per sample")
        plt.xlabel("Samples")
        plt.ylabel("Total no. of numts")
        plt.ylim(0,370)
        # plt.yticks(range(0,370,10))
        plt.savefig("Num_of_total_numts_every_sample.png")
        plt.close()
        return all_files

    def natural_sort(self,l):
        # a = [i.spl for i in l]
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alpha_num = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key= alpha_num)

    def get_df_from_pickle(self):
        all_unique_numts = defaultdict(list)
        for idx, subdir in enumerate(self.flank_df):
            s_dir1 = os.listdir(os.path.join(hardcoded_path.flank_df,subdir))
            s = self.natural_sort(s_dir1)
            # print(s)
            # print(s_dir1[0].split(".")[0][3:])
            # print(s)
            # all_unique_numts[subdir_no] = unique_numts_file

            for file in s:
            #     # print(file)
                read_pickle_file = pd.read_pickle(os.path.join(hardcoded_path.flank_df,subdir,file))
                # print(read_pickle_file)
            # print("blahhhh")
                
            # break
                unique_numts = read_pickle_file['Numt'].unique()
                # print(len(unique_numts))
                all_unique_numts[subdir].append(len(unique_numts))
                # for i in unique_numts:
                # #     #How many numts with 1 supporting read were identified
                #     if i in self.numts_1read[idx]:
                # #         # print("blahblah")
                #         print(i)
                # print("blah")
                # print(unique_numts)
                # break

        # print(all_unique_numts)
        # idx = 0
        # for k, v in all_unique_numts.items():
        #     plt.bar(range(len(v)), v)
        #     for index,data in enumerate(v):
        #         plt.text(x=index , y =data , s=f"{data}" , fontdict=dict(fontsize=10), ha = 'center')
        #     plt.title(f"No. of numts per kmer for {k} out of {self.all_files_total_nmts[idx]} actual numts")
        #     plt.xlabel("ksize")
        #     plt.ylabel("No. of numts")
        #     plt.ylim(0,360)
        #     plt.xticks(range(len(v)), range(17,36,2))
        #     # plt.show()
        #     plt.savefig(f"No_numts_allkmers_persample_{k}.png")
        #     plt.close()
        #     idx += 1
        return all_unique_numts

    def num_reads_mt_origin(self):
        all_files = []
        for idx, name in enumerate(self.mt_reads):
            num_mt_reads = [] 
            sam_file = pysam.AlignmentFile(os.path.join(hardcoded_path.Mt_reads_dir,name),'rb')
            # print(name)
            for aln in sam_file:
                num_mt_reads.append(aln.query_name)

            # print(name.split(".")[0],"No. of reads that are Mt in origin" ,len(num_mt_reads))
            all_files.append((name.split(".")[0],len(num_mt_reads)))
        # print(all_files)
        return all_files

    def num_reads_candidate_nmts(self):
        all_files = []
        for idx, name in enumerate(self.numt_reads):
            num_numt_reads = []
            sam = pysam.AlignmentFile(os.path.join(hardcoded_path.numt_reads,name),'rb')
            for aln in sam:
                num_numt_reads.append(aln.query_name)

            # print(name.split(".")[0],"No. of reads that are Mt in origin" ,len(num_numt_reads))
            all_files.append((name.split(".")[0],len(num_numt_reads)))
        # print(all_files)
        return all_files

    def get_result_df(self):
        ovlp_data = self.num_reads_ovlp_MT_nmts()
        df = pd.DataFrame(ovlp_data, columns =['File', 'No. Reads overlapping Numts/MT'])
        mt_data = self.num_reads_mt_origin()
        df1 = pd.DataFrame(mt_data, columns =['File', 'No. reads Mt in origin'])
        numt_data = self.num_reads_candidate_nmts()
        df2 = pd.DataFrame(numt_data, columns =['File', 'No. candidate numt reads'])
        all_df = df.merge(df1,on='File').merge(df2,on='File')
        print(all_df)
        return all_df

    def remove_dict_key(self, my_dict, num_files):
        y = dict()
        for k,v in my_dict.items():
            if v == num_files:
                y[k] = v
        return y


    def get_sample_specific_or_fixed(self):
        # self.unique_nm
        specific_numts_df = defaultdict(dict)
        fixed_numts_df = defaultdict(dict)
        for all_files in self.flank_df:
            print(all_files)
            # unique_nmt = []
            read_pickle_file = pd.read_pickle(os.path.join(hardcoded_path.flank_df,all_files))
            unq_numts = list(read_pickle_file['Numt'].unique())
            val_not_in_all_df = defaultdict(lambda: 0)
            val_in_all_df = defaultdict(lambda: 0)
            for remaining_files in self.flank_df:
                if all_files != remaining_files:
                    read_pickle_again = pd.read_pickle(os.path.join(hardcoded_path.flank_df,remaining_files))
                    # for rows in read_pickle_again.itertuples():
                    unq_numts_again = list(read_pickle_again['Numt'])
                    for val in unq_numts:
                    # print("red pkl: ", read_pickle_again['Numt'])
                        if val not in unq_numts_again:
                            val_not_in_all_df[val] += 1
                        elif val in unq_numts_again:
                            val_in_all_df[val] += 1
            only_specific = self.remove_dict_key(val_not_in_all_df, len(self.flank_df) - 1)
            specific_numts_df[all_files] = only_specific
            # print(val_not_in_all_df)
            
            # print(val_in_all_df)
        only_fixed = self.remove_dict_key(val_in_all_df, len(self.flank_df) - 1)
        fixed_numts = [k for k in only_fixed.keys()]
            # fixed_numts_df[all_files] = len(only_fixed)
        # print("Specific_numts",specific_numts_df)
        # print("fixed numts", fixed_numts_df)
                # break
            # break
        return specific_numts_df, fixed_numts

    def get_all_numts_for_all_samples(self):
        c = 0
        # file1_ = []
        all_ref_numts = []
        for line in open(hardcoded_path.bed_dir_without_chrM):
                # pd.DataFrame(ovlp_data, columns =['File', 'No. Reads overlapping Numts/MT'])
                fields = line.rstrip().split("\t")
                all_ref_numts.append(fields[3])
        df = pd.DataFrame(all_ref_numts, columns =['Numt'])
        all_sample_names = []
        # print(df)
        sample = defaultdict(list)
        for files in self.flank_df:
            read_pickle_file = pd.read_pickle(os.path.join(hardcoded_path.flank_df,files))
            unq_numts = list(read_pickle_file['Numt'].unique())
            for nmt in all_ref_numts:
                if nmt in unq_numts:
                    sample[files].append(1)
                    # print(fields[3])
                    #append 1 if fields[3] in unq_numts
                else:
                    sample[files].append(0)
                   
        for k,v in sample.items():
            s = k.split(".")[0]
            all_sample_names.append(s)
            df[f"{s}"] = v
        # print(df)
        # print(c)
        # pd.set_option('display.max_rows', None) 
        # pd.set_option('display.max_colwidth', 16)
        # print(all_sample_names)
        counts_df = df.groupby(all_sample_names).size().reset_index().sort_values([0], ascending= False).rename(columns={0:'count'})
        # print(counts_df['count'])
        # df.set_index('Numt', inplace= True)
        df['Presence_Freq'] = df.sum(axis=1, numeric_only = True) / len(all_sample_names)
        df['Absence_Freq'] = 1 - df['Presence_Freq']
        freq_list = df['Presence_Freq'].to_list()
        plt.hist(freq_list)
        # hist = df['Frequency'].hist()
        plt.savefig("Sample_freq.png")
        # print(hist)
        # df['Count'] = df.sum(axis=1, numeric_only = True) 
        
        # print(df)
        pd.set_option('display.max_rows', None) 
        # pd.set_option('display.max_colwidth', None)
        # print(counts_df)
        return df, counts_df

# result = get_result()
# print(result.get_all_numts_for_all_samples())

# print(result.get_sample_specific_or_fixed())
# result.get_result_df()
# result.num_reads_ovlp_MT_nmts()
# result.get_df_from_pickle()
# result.num_reads_mt_origin()
# result.num_reads_candidate_nmts()