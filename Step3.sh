#!/bin/sh
OUTFILE=""
INFILE=""
reference=""
for i in "$@"; do
    case $i in
        -o=*|--output=*)
          OUTFILE="${i#*=}"
          shift
          ;;
        -i=*|--input=*)
          INFILE="${i#*=}"
          shift
          ;;
        -ref=*|--reference=*)
          reference="${i#*=}"
          shift
          ;;
        *)
        ;;
    esac
done

# samtools index original_files/outfile/ERR1019039.MT.csort.alignment.numt_tag.bam
samtools index "${INFILE}"
##So you give the reference sequence (human genome .fasta file) and the sample bamfile to bcftools mpileup and its outputs a .bcf file.
# bcftools mpileup -Ou -r chrM -f original_files/Reference_Genome_GRCh38.p13/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa original_files/outfile/"${INFILE}"  > original_files/Temp/Temp_chrM.bcf
bcftools mpileup -Ou -r chrM -f "${reference}" "${INFILE}"  > Temp/Temp_chrM.bcf

##Pipe the output of bcftools mpileup to bcftools call to call the variants stored in bcf file and give a .vcf file
bcftools call -v -c -Oz --ploidy 1 -o Temp/Temp_chrM.vcf.gz Temp/Temp_chrM.bcf

##For consensus sequence I need to get the index of the vcf file
bcftools index Temp/Temp_chrM.vcf.gz

samtools faidx "${reference}" chrM | bcftools consensus Temp/Temp_chrM.vcf.gz > Temp/consensus_chrM_new.fa

##Sort first
# samtools sort -o Temp/ERR1019039.numt_tag.sort.bam -n original_files/outfile/ERR1019039.MT.csort.alignment.numt_tag.bam
samtools sort -o Temp/Temp.numt_tag.sort.bam -n "${INFILE}"
# samtools sort -o Temp/Temp.numt_tag.sort.bam "${INFILE}"

##Get paired fq files from the numt_tag_new.bam
samtools fastq -1 Temp/numt_tag_R1.fq -2 Temp/numt_tag_R2.fq Temp/Temp.numt_tag.sort.bam
# samtools fastq  ERR1019039.MT.csort.alignment.numt_tag_new.bam

##Index the ref.fa file 
bwa index Temp/consensus_chrM_new.fa

##Use that .fq file along with Consensus Mt as Reference and get the .bam file as output. 
#bwa mem -o original_files/Step3_output/ERR1019039.MT.Genome.numt_tag.consensus_Mt.bam Temp/consensus_chrM_new.fa Temp/numt_tag_R1.fq Temp/numt_tag_R2.fq
bwa mem -o "${OUTFILE}" Temp/consensus_chrM_new.fa Temp/numt_tag_R1.fq Temp/numt_tag_R2.fq

rm -r Temp
