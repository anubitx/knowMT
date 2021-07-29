import os

current_dir = os.getcwd()
##original_files and bed_file are the folders that I have created
orignal = os.path.join(current_dir, 'original_files')
bed_dir = os.path.join(current_dir, 'original_files', 'bed_file_with_chrM')
out_dir = os.path.join(current_dir, 'original_files', 'outfile')
# bed_dir = os.path.join(bed_dir, os.listdir(bed_dir)[0])
# bed_dir_without_chrM =  os.path.join(current_dir, 'original_files', 'bed_file')
# bed_dir_without_chrM = os.path.join(bed_dir_without_chrM, os.listdir(bed_dir_without_chrM)[0])
bam_dir = os.path.join(current_dir, 'original_files', 'bam_files')
# temp_dir = os.path.join(current_dir, 'original_files', 'Temp')
# step3_output_dir = os.path.join(current_dir, 'original_files', 'Step3_output')
# stepf_outfile = os.path.join(current_dir, 'original_files', 'Stepf_outfile')
# Mt_reads_dir = os.path.join(current_dir, 'original_files', 'Mitochondrial_originating_reads')
# numt_reads = os.path.join(current_dir, 'original_files', 'Only_numt_reads')
# chrM_str_path = os.path.join(current_dir, 'chrM_string.txt')
# kmer_fastq_files = os.path.join(current_dir, 'original_files', 'kmer_fastq_files')
# Mt_input_bam = os.path.join(current_dir, 'original_files','Mt_input_bam')
# Minimap_output = os.path.join(current_dir, 'original_files','Minimap_output')

# make_dir_list = [orignal,out_dir,temp_dir,step3_output_dir,stepf_outfile,Mt_reads_dir, numt_reads,kmer_fastq_files, Mt_input_bam,Minimap_output]
make_dir_list = [orignal,out_dir,bed_dir,bam_dir]
def make_dirs():
    for i in make_dir_list:
        if not os.path.exists(i):
            os.mkdir(i)
            
# Spades_output_ERR1019039
# print(Spades_output)

# print(Mt_input_bam)
# print(current_dir)
# print(bed_dir)
# print(bam_dir)
# print(kmer_fastq_files)