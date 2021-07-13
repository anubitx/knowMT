import pysam
import os
import simplesam
from intervaltree_bio import GenomeIntervalTree
import hardcoded_path
import subprocess
from tqdm import tqdm

## Checking for the dna fragment instead of individual reads for overlap

##Conducted a test to cross check if the operation that consumes the query is accounted for or not. 
##It is accounted for, hence we use simplesam package.
# y = simplesam.Sam(cigar = '5M5I5N5=')
# print(len(y))

# gr_read = pr.PyRanges(chromosomes = [read.reference_name], starts= [read.reference_start], ends = [read.reference_end])
# gr_mate = pr.PyRanges(chromosomes = [read.next_reference_name], starts= [read.next_reference_start], ends = [read.next_reference_start + len(x)])
    
class numt_chrM_overlapping_reads():
    def __init__(self, bed_file_path:str, bam_file_path:str ) -> None:
        
        self.ref_numt = GenomeIntervalTree()
        self.bed = bed_file_path
        self.bam_folder = bam_file_path
        self.read_bed()
        self.bam_file_name = self.get_bam()
        self.outfiles = self.parse_bam_name_for_outfile()
        self.step3_out = self.parse_bam_name_for_step3_outfile()
        self.generate_numt_tag()
        # self.step3_for_all_files()
        # print(self.outfiles)

    def parse_bam_name_for_outfile(self):
        out_name = "numt_tag"
        outfile_list = []
        for bam_name in self.bam_file_name:
            s = bam_name.split('.')
            s.insert(-1, out_name)
            outfile_list.append('.'.join(s))
        return outfile_list
    
    def parse_bam_name_for_step3_outfile(self):
        out_name = ".Genome.numt_tag.consensus_Mt.bam"
        outfile_list = []
        for bam_name in self.outfiles:
            s = bam_name.split('.')
            s = f"{s[0]}.{s[1]}{out_name}"
            outfile_list.append(s)
        return outfile_list

    def read_bed(self):
        for line in open(self.bed):
            fields = line.rstrip().split("\t")
            self.ref_numt[fields[0]].addi(int(fields[1]), int(fields[2]), (fields[3], int(fields[2]) - int(fields[1])))

    def get_bam(self):
        return os.listdir(self.bam_folder)

    def numt_Overlap(self, chr:str, start, end):
        overlap = 0
        overlapid = "Mapped"
        overlapcontained = 0
        overlap_count = 0
        nlen = -1
        refoverlap = self.ref_numt[chr].search(start,end)
        if len(refoverlap) > 0:
            overlap = 1
            overlap_count = len(refoverlap)
            refoverlap = refoverlap.pop()
            overlapid = refoverlap[2][0]
            nlen = refoverlap[2][1]
            if refoverlap[0] <= start and refoverlap[1] >= end:
                overlapcontained = 1
            elif refoverlap[0] >= start and refoverlap[1] <= end:
                overlapcontained = 2      
        return str(overlap), str(overlapcontained), overlapid, str(nlen) , str(overlap_count)

    def generate_numt_tag(self):
        for idx,i in enumerate(tqdm(self.bam_file_name)):
            samfile = pysam.AlignmentFile(os.path.join(self.bam_folder,i), 'rb')
            outfile = pysam.AlignmentFile(os.path.join(hardcoded_path.out_dir, self.outfiles[idx]), 'wb', header = samfile.header)
            reads_retained = []

            for read in samfile:
                if (not(read.is_unmapped)):
                    self_info = self.numt_Overlap(read.reference_name, read.reference_start, read.reference_end)
                else:
                    self_info = ("0", "0", "Unmapped", "-1" , "0")
                if (not( read.mate_is_unmapped)):
                    x = simplesam.Sam(cigar = read.get_tag("MC"))
                    mate_info = self.numt_Overlap(read.next_reference_name, read.next_reference_start, read.next_reference_start + len(x))
                else:
                    mate_info = ("0", "0", "Unmapped", "-1" , "0")
                if self_info[0] == "1" or mate_info[0] == "1":
                    reads_retained.append(read.query_name)
                    read.set_tag("Ni", ','.join(self_info + mate_info), value_type = "Z", replace = False)
                    outfile.write(read)
            print("Total number of reads that overlap numts/MT:",len(reads_retained))
            outfile.close()  

    def step3_for_all_files(self):
        # p = "/mnt/c/Users/simar/JCSMR_numtS_Project/ANU_Numts/nmt/Step3.sh"
        cwd = os.getcwd()
        p = os.path.join(cwd, "Step3.sh")
        for idx, file_name in enumerate(self.outfiles):
            in_path = os.path.join(hardcoded_path.out_dir, file_name)
            out_path = os.path.join(hardcoded_path.step3_output_dir, self.step3_out[idx])
            print(in_path,out_path)
            cmd = f" -o={out_path} -i={in_path}"
            cmd = p+cmd
            print(cmd)
            s = subprocess.Popen([cmd], shell = True)
            out, err = s.communicate()



# numt_chrM_overlapping_reads.read_bed(bed, bam)
