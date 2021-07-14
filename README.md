# nmt
Discovery and annotation of polymorphic nuMTs in the human genome.  

Dependencies  
Python 3.6  
pip install pysam  
pip install intervaltree_bio  
pip install simplesam  
pip install tqdm  

1. From the git repository Download the following files:  
Numt_Tag.py  
hardcoded_path.py  
main.py  
Ref_numt_grouped_blocks_final_with_chrM.bed  

2. Install all the above mentioned dependencies

3. Run the following command:  
python3 main.py -bam-folderpath <path/to/bam/your/folder> -bed-filename <path/to/the/downloaded/bed/file/Ref_numt_grouped_blocks_final_with_chrM.bed>  

Note:  
-bam-folderpath <path/to/bam/folder> is the path to the folder that contains input bam files.  
-bed-filename <path/to/bed/filename> is the path to the input bed file.  

4. The output of this command is stored in a folder named outfile.   