# knowMT
## KnowMT: a software to detect Nuclear insertions of Mitochondrial Sequences (NumtS) in humans using whole genome sequence data  
I have developed knowMT software to detect NumtS in the human whole genome sequence data. KnowMT successfully removes mitochondrial noise that confounds the detection of NumtS in the whole genome sequence data by (a) generating individual specific Mt sequence to filter Mt reads and (b) mask Mt similar sequences from Numts. Using this method, the flanking regions of the NumtS are enriched and detected to estimate the insertion points of NumtS in the genome.   

Dependencies  
Python 3.6  
pip install pysam  
pip install PyIntervalTree  
pip install simplesam  
pip install tqdm  

1. From the git repository Download the following files:  
Numt_Tag.py  
var_ovlp.py
alignment_approach.py
cluster.py
main.py  
Ref_numt_grouped_blocks_final_with_chrM.bed  

2. Install all the above mentioned dependencies

3. Parameters:  
-input_file default None,
-numt_bed default None,
-out_dir default None,
-reference default None,
-threads default 1,
-ksize default 35,
-path_to_spades None

4. Run the following command:  
python3 main.py -input_file <path/to/bam/your/file> -numt_bed <path/to/the/downloaded/bed/file/Ref_numt_grouped_blocks_final_with_chrM.bed> -out_dir <path/to/output/directory> -reference <path/to/reference.fasta/file> -threads <integer> -ksize <kmer/size/to/mask/reads> -path_to_spades <path/to/the/spades.py>  

5. The output of this command is stored in a folder named OUTPUT.   
