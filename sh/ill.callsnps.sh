#!/bin/bash
# Map reads to phage genome
numcpu=$1
for sample in `ls reads.ill | grep R1 | sed 's/\.R1\.fastq//g'`
do
	bowtie2 -p $numcpu --maxins 8326 -x SP055-rpoZ-cMyc-Cry1Ac1-d123 -1 reads.ill/${sample}_R1.fastq -2 reads.ill/${sample}_R2.fastq -S map.ill/$sample.sam
done

# Convert to BAM
for sample in `ls map.ill | sed 's/\.sam//g'`
do
	samtools view -bS map.ill/$sample.sam > map.ill/$sample.bam
done

# Call variants
cd map.ill
cmd=bamaddrg
for sample in `ls | sed 's/\..am//g' | uniq`
do
	cmd="$cmd -b $sample.bam -s $sample -r $sample"
done
`$cmd > allbam.bam`
samtools sort allbam.bam allbam.sorted
samtools rmdup allbam.sorted.bam allbam.sorted.rmdup.bam
samtools index allbam.sorted.rmdup.bam
outfile=freebayes.`date +%Y%m%d%H%M%S`.vcf
freebayes --use-best-n-alleles 1 --pooled-continuous --use-reference-allele --theta 500000000 --min-alternate-fraction 0.01 --ploidy 1 --region SP055-rpoZ-cMyc-Cry1Ac1-d123:2833-4971 -f ../SP055-rpoZ-cMyc-Cry1Ac1-d123.fasta allbam.sorted.rmdup.bam > $outfile

# Prepare table for analysis in R
outpfx=`echo $outfile | sed 's/.vcf//g'`
grep -v '^##' $outfile | sed 's/#//g' > $outpfx.data.tsv
