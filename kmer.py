def read_input_file(filename):
    with open(filename, 'r') as f:
        dna = f.read()
    return dna

# Mt_string = read_input_file("chrM_string.txt")
# print(Mt_string)

def reverse_complement(dna_input):
    """return the REVERSE complement of the dna string provided """
    reverse = dna_input[::-1].strip()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = [complement[base] for base in reverse] 
    return ''.join(bases)

# Mt_Rc = reverse_complement(Mt_string)
# print(Mt_Rc)

##Function that takes a sequence and ksize and returns kmers of ksize from the given sequence and its reverse complement. 
def build_kmers(sequence, ksize):
    kmers = []
    Mt_rc = reverse_complement(sequence)
    # print(type(Mt_rc))
    n_kmers = len(sequence) - ksize + 1
    n_rc_kmers = len(Mt_rc) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    for i in range(n_rc_kmers):
        kmer = Mt_rc[i:i + ksize]
        kmers.append(kmer)

    return kmers


import pysam

bamfile = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_tag_new.bam', 'rb')
# fasta_file = open('ERR1019039.numt_tag_maskedN.fasta', 'w')

Mt_string = read_input_file("chrM_string.txt")
Mt_kmers = set(build_kmers(Mt_string, 31))
# print(Mt_kmers)
# counter = 0

# KSIZE = 31
def check_kmers(sequence, ksize):
    modseq = ""
    #n_kmers = len(sequence) - ksize + 1
    for i in range(len(sequence) - ksize + 1):
        kmer = sequence[i:i + ksize]
        if kmer in Mt_kmers:
            kmer = 'N' * ksize
        if i > 0:
            modseq = modseq[:i] + kmer  
        else:
            modseq = kmer
    return modseq

for aln in bamfile:
    # print(aln.query_name)

    if aln.is_read1 == True:
        # print('>' + aln.query_name + '/1' + '\n' + check_kmers(aln.query_sequence, 31) + '\n')
        fasta_file.write('>' + aln.query_name + '/1' + '\n' + check_kmers(aln.query_sequence, 31) + '\n')
    elif aln.is_read2 == True:
        # print('>' + aln.query_name + '/2' + '\n' + check_kmers(aln.query_sequence, 31) + '\n')
        fasta_file.write('>' + aln.query_name + '/2' + '\n' + check_kmers(aln.query_sequence, 31) + '\n')
   
   
# fasta_file.close()
   
   
   
   
   
# file1 = open('chrM.fasta', 'r')

# lines = file1.readlines()
# file1.close()
# new_file = open("chrM.txt", "w")
# for line in lines:
#     if line.strip("\n") != ">MT":
#         new_file.write(line)
# new_file.close()

# file1 = open('chrM.txt', 'r')
# new_file = open('chrM_kmers.txt', 'w')
# string_without_line_breaks = ""
# for line in file1:
#     stripped_line = line.rstrip()
#     string_without_line_breaks += stripped_line
# 
#    
    # print(read_kmers)
    ## Intersection of both the kmer sets, to find out the 31-mers that are common in the mitochondrial sequence and the query_sequence
    #result = Mt_kmers.intersection(read_kmers)
    #print(result)
    # print(aln , aln.query_sequence)

    ##Initialise lists within the loop, else the minimum and maximum of all the alignments will be stored in one list.
    ##And we will get only one minimum and maximum for all the alignments. 
    ##To avoid that make 2 temporary lists within the loop. 
#     temp_min_list = []
#     temp_max_list = []

#     ##For each 31-mer in the result. (ie. all the 31-mers that are common in mitochondrial sequence and the query sequence)
#     for i in result:
#         ## Replace the string i (every common 31-mer) with len(i) times N. 
#         ## Else string of len 31 will be replaced with 1 N.
#         a = aln.query_sequence.replace(i, 'Z' * len(i))
#         # print(a)
#         # print(a)
#         ## Get the index of first occurence of 'N' for all 31-mers in the result
#         ## If there is an exisiting N in the query sequence then it gives me that index as the minimum 
#         ## which is not what I want
#         temp_min_n = a.index('Z')
#         ## Append these values to a temporary list and 
#         temp_min_list.append(temp_min_n)
#         temp_max_n = a.rfind('Z')
#         temp_max_list.append(temp_max_n)
#         # min_n = len(a)
#         # print(temp_min_n, temp_max_n)
    
#     ##If the list is not empty then get the minimum and maximum values.
#     if len(temp_min_list) != 0:
#         min_val = min(temp_min_list)
#     if len(temp_max_list) != 0:
#         max_val = max(temp_max_list)
#         # print(min_val, max_val)
#     ##If the result is an empty set just return the query sequence
#     if len(result) == 0:
#         string = aln.query_sequence
#         # print( string)
#         # print(string, aln)
#     ##If the result is not a non empty set then replace the strings matching to all elements of result with N
#     if len(result) > 0:
#         string = aln.query_sequence
#         # print(string, len(string))
#         string = string[:min_val] + 'Z' * (max_val - min_val + 1) + string[max_val+1:]
#         # print(string.replace('Z', 'N'))
#         aln.query_sequence = string.replace('Z', 'N')
#         print( aln.query_sequence)
#         # print(final)
#         # print(string, aln)
# # print(string)
       
    # counter += 1
    # if counter == 2200:
        # break
    
