#!/bin/bash
topDir=$(pwd)

####STAR paramaters
STAR_runThreadN=1
STAR_outFileNamePrefix=${topDir}"/result"
STAR_readFilesCommand=""

while getopts "g:t:o:za:b:n:T:f:c:m:" opt; do
	case $opt in
	g) STAR_genomeDir=$OPTARG ;;
	t) STAR_runThreadN=$OPTARG ;;
	o) STAR_outFileNamePrefix=$OPTARG ;;
	z) STAR_readFilesCommand="zcat" ;;
	a) STAR_readFilesIn_1=$OPTARG ;;
	b) STAR_readFilesIn_2=$OPTARG ;;
	n) STAR_outFilterMatchNmin=$OPTARG ;;
	T) TE_annotation_bedfile=$OPTARG ;;
	f) genome_fa_file=$OPTARG ;;
	c) genome_chr_file=$OPTARG ;;
	m) mappability_bedGraph_file=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

#### map raw reads
mkdir ${STAR_outFileNamePrefix}
STAR   --runThreadN $STAR_runThreadN --genomeDir ${STAR_genomeDir}   --readFilesIn ${STAR_readFilesIn_1} ${STAR_readFilesIn_2}  --readFilesCommand $STAR_readFilesCommand  --outReadsUnmapped Fastx   --outFilterMultimapNmax 1   --outFilterMatchNmin $STAR_outFilterMatchNmin   --outFilterMismatchNmax 2   --outFilterMismatchNoverLmax 0.05 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${STAR_outFileNamePrefix}
bamToBed -i ${STAR_outFileNamePrefix}/Aligned.sortedByCoord.out.bam -split >${STAR_outFileNamePrefix}/Aligned.out.bed

#### sort split alignments
python ${topDir}/scripts/split_sort.py ${STAR_outFileNamePrefix} >${STAR_outFileNamePrefix}/line1.txt
sortBed -i ${STAR_outFileNamePrefix}/line.txt >${STAR_outFileNamePrefix}/line.sorted.txt
mergeBed -i ${STAR_outFileNamePrefix}/line.sorted.txt -c 1 -o count >${STAR_outFileNamePrefix}/line.merge.txt
mergeBed -i ${STAR_outFileNamePrefix}/line1.txt -c 1 -o count >${STAR_outFileNamePrefix}/line1.merge.txt
cat ${STAR_outFileNamePrefix}/line.merge.txt ${STAR_outFileNamePrefix}/line1.merge.txt >${STAR_outFileNamePrefix}/line_split.txt
sortBed -i ${STAR_outFileNamePrefix}/line_split.txt >${STAR_outFileNamePrefix}/line_split.sorted.txt
mergeBed -i ${STAR_outFileNamePrefix}/line_split.sorted.txt -c 4 -o sum >${STAR_outFileNamePrefix}/Aligned.out.merge.bed
rm ${STAR_outFileNamePrefix}line.txt ${STAR_outFileNamePrefix}/line1.merge.txt ${STAR_outFileNamePrefix}/line1.txt ${STAR_outFileNamePrefix}/line.sorted.txt ${STAR_outFileNamePrefix}/line.merge.txt ${STAR_outFileNamePrefix}/line_split.txt ${STAR_outFileNamePrefix}/line_split.sorted.txt

####uniq 
intersectBed -a ${STAR_outFileNamePrefix}/Aligned.out.merge.bed -b ${TE_annotation_bedfile} -wa -wb > ${STAR_outFileNamePrefix}/RMSK_uniq_reads_overlaped_filtered_merge.bed
python ${topDir}/scripts/TE_transcript_grow.py ${STAR_outFileNamePrefix}/RMSK_uniq_reads_overlaped_filtered_merge.bed >${STAR_outFileNamePrefix}/RMSK_transcript.txt
sortBed -i ${STAR_outFileNamePrefix}/RMSK_transcript.txt >${STAR_outFileNamePrefix}/RMSK_transcript.sorted.txt
mergeBed -i ${STAR_outFileNamePrefix}/RMSK_transcript.sorted.txt  -c 1,4,6,5 -o count,distinct,sum,distinct >${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.txt
rm ${STAR_outFileNamePrefix}/RMSK_uniq_reads_overlaped_filtered_merge.bed ${STAR_outFileNamePrefix}/RMSK_transcript.txt ${STAR_outFileNamePrefix}/RMSK_transcript.sorted.txt

####get TE sequence that overlaped with uniq alignments
awk 'BEGIN {id="1"}{if ($7~/,/)$7="+"} {print $1"\t"$2"\t"$3"\t"id"\t"$5"\t"$7}{id=id+1}' ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.txt >${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.bed
rm ${STAR_outFileNamePrefix}/mm9_rmsk.fasta
mkdir ${STAR_outFileNamePrefix}/fastafrombed
for i in $(cat ${genome_chr_file})
do    
	grep -w $i ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.bed >${STAR_outFileNamePrefix}/fastafrombed/$i.uniq.bed    
	fastaFromBed -fi ${genome_fa_file}/$i.fa -bed ${STAR_outFileNamePrefix}/fastafrombed/$i.uniq.bed -fo ${STAR_outFileNamePrefix}/fastafrombed/$i.fasta -s -name    
	cat ${STAR_outFileNamePrefix}/fastafrombed/$i.fasta >>${STAR_outFileNamePrefix}/mm9_rmsk.fasta
done

####getfasta correction for bedtools bug
python ${topDir}/scripts/getfasta_correction.py ${STAR_outFileNamePrefix}/mm9_rmsk.fasta >${STAR_outFileNamePrefix}/mm9_rmsk.fa
mv ${STAR_outFileNamePrefix}/mm9_rmsk.fa ${STAR_outFileNamePrefix}/mm9_rmsk.fasta

####build TE sequence index
mkdir ${STAR_outFileNamePrefix}/STARIndex
STAR --runThreadN $STAR_runThreadN --runMode genomeGenerate --genomeDir ${STAR_outFileNamePrefix}/STARIndex --genomeFastaFiles ${STAR_outFileNamePrefix}/mm9_rmsk.fasta --genomeChrBinNbits 10  #--genomeChrBinNbits set to avoid too much memoery used

####get mappability file
intersectBed -a ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.bed -b ${mappability_bedGraph_file} -wa -wb >${STAR_outFileNamePrefix}/RMSK_transcript_mappability_intersectbed.bed
python ${topDir}/scripts/process_mappability.py ${STAR_outFileNamePrefix}/RMSK_transcript_mappability_intersectbed.bed >${STAR_outFileNamePrefix}/mappability_score_of_TE.txt
rm ${STAR_outFileNamePrefix}/RMSK_transcript_mappability_intersectbed.bed

####map unmapped pair-end reads separetly
mkdir ${STAR_outFileNamePrefix}/Unmapped_mate1/
mkdir ${STAR_outFileNamePrefix}/Unmapped_mate2/
STAR --runThreadN $STAR_runThreadN --alignEndsType EndToEnd --genomeDir ${STAR_outFileNamePrefix}/STARIndex/ --readFilesIn ${STAR_outFileNamePrefix}/Unmapped.out.mate1 --outFilterMatchNminOverLread 1 --scoreDelOpen -100 --scoreDelBase -100 --scoreInsOpen -100 --scoreInsBase -100 --scoreGapNoncan -100 --scoreGap -100 --outFileNamePrefix ${STAR_outFileNamePrefix}/Unmapped_mate1/
STAR --runThreadN $STAR_runThreadN --alignEndsType EndToEnd --genomeDir ${STAR_outFileNamePrefix}/STARIndex/ --readFilesIn ${STAR_outFileNamePrefix}/Unmapped.out.mate2 --outFilterMatchNminOverLread 1 --scoreDelOpen -100 --scoreDelBase -100 --scoreInsOpen -100 --scoreInsBase -100 --scoreGapNoncan -100 --scoreGap -100 --outFileNamePrefix ${STAR_outFileNamePrefix}/Unmapped_mate2/ 

####calculate TE expression
awk 'BEGIN {id="1"}{print id"\t"$6}{id=id+1}' ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.txt >${STAR_outFileNamePrefix}/RMSK_transcript_uniq.cntTable
python ${topDir}/scripts/ZX.py ${STAR_outFileNamePrefix}/Unmapped_mate1/Aligned.out.sam ${STAR_outFileNamePrefix}/mappability_score_of_TE.txt ${STAR_outFileNamePrefix}/RMSK_transcript_uniq.cntTable > ${STAR_outFileNamePrefix}/Unmapped_mate1/RMSK_transcript.cntTable
python ${topDir}/scripts/ZX.py ${STAR_outFileNamePrefix}/Unmapped_mate2/Aligned.out.sam ${STAR_outFileNamePrefix}/mappability_score_of_TE.txt ${STAR_outFileNamePrefix}/RMSK_transcript_uniq.cntTable > ${STAR_outFileNamePrefix}/Unmapped_mate2/RMSK_transcript.cntTable

####Integrate results		
python ${topDir}/scripts/combined_uniq_multi_count.py ${STAR_outFileNamePrefix}/RMSK_transcript_uniq.cntTable ${STAR_outFileNamePrefix}/Unmapped_mate1/RMSK_transcript.cntTable >${STAR_outFileNamePrefix}/RMSK_transcript_uniq_mate1.cntTable
python ${topDir}/scripts/combined_uniq_multi_count.py ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_mate1.cntTable ${STAR_outFileNamePrefix}/Unmapped_mate2/RMSK_transcript.cntTable >${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.cntTable
sort -n ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.cntTable >${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.sorted.cntTable
cut -f2 ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.sorted.cntTable >${STAR_outFileNamePrefix}/b
paste ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.bed ${STAR_outFileNamePrefix}/b >${STAR_outFileNamePrefix}/RMSK_transcript.cntTable
rm ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_mate1.cntTable ${STAR_outFileNamePrefix}/b ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.sorted.cntTable ${STAR_outFileNamePrefix}/RMSK_transcript_uniq_multi.cntTable ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.bed ${STAR_outFileNamePrefix}/RMSK_transcript.aggrated.txt ${STAR_outFileNamePrefix}/RMSK_transcript_uniq.cntTable
