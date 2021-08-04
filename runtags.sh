#!/bin/bash
#PBS -j oe
#PBS -N filterbams
#PBS -o /g/data/te53/numts/logs/
#PBS -l ncpus=48,walltime=48:00:00,mem=190GB

module load python3packages parallel

outputdir=/g/data/te53/numts/filteredbams/
ref=/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa
numtbed=/g/data/te53/numts/nmt/Ref_numt_grouped_blocks_final_with_chrM.bed
parallellog=/g/data/te53/numts/filteredbams/numttags_parallel.H$head.T$tail.log

ls /g/data/te53/phase2_20200312/SAMEA*/realignedcram/*.recal.realn.cram |\
head -n $head | tail -n $tail |\
parallel --resume-failed --results /g/data/te53/hrp561/covseq/alignments/cmdlog/ --joblog $parallellog -j 12 python3 /g/data/te53/numts/nmt/main.py -inputbam {} -numtbed $numtbed -outputdir $outputdir -threads 4 -reference $ref


