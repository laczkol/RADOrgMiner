#!/bin/bash
set -euo errexit -euo pipefail -euo nounset

echo "This is RADOrgMiner (RADseq Organellar DNA Miner and Genotyper) v0.8"

echo "Run started at $(date)"

cdir=$(pwd)
indir=$(pwd)
outdir=$(pwd)
ref=" "
mask="no"
np=1
map="yes"
popmap=" "
call="no"
type="SE"
bwa_k=19
bwa_A=1
bwa_B=4
bwa_O=6
min_bed_cov=3
sub_prop=1.0
ploidy=1
minqual=20
mismatchqual=15
mapq=30
fb_min_cov=5
nbest=5
altfrac=0.4
vcfmiss=0.8
mbs=0
mrc=1
minlen=100

while [[ "$#" -gt 0 ]];
	do
	case $1 in
		-s|--samples-directory)
			indir="$2"
			shift
			;;
		-o|--output-directory)
			outdir="$2" #output directory, USE ABSOLUTE PATH
			shift
			;;
		-r|--reference-genome) #reference genome in fasta, MANDATORY
			ref=$2
			shift
			;;
		-mask|--mask-reference) #specify yes to mask genome
			mask=$2
			shift
			;;
		-np|--number-of-processors) #number of processors
			np=$2
			shift
			;;
		-align|--align-reads) #yes or no
			map=$2
			shift
			;;
		-bwa_k|--min-seed-length) #minimum seed length of bwa
			bwa_k=$2
			shift
			;;
		-bwa_A|--matching-score) #matching score for bwa
			bwa_A=$2
			shift
			;;
		-bwa_B|--mismatch-penalty) #mismatch-penalty for bwa
			bwa_B=$2
			shift
			;;
		-bwa_O|--gap-open-penalty) #gap-open penalty for bwa
			bwa_O=$2
			shift
			;;
		-popmap|--population-map) #population map MANDATORY
			popmap=$2
			shift
			;;
		-type|--sequencing-type) #SE or PE
			type=$2
			shift
			;;
		-call|--call-haplotypes) #yes or no
			call=$2
			shift
			;;
		-mbc|--min-bed-coverage)
			min_bed_cov=$2
			shift
			;;
		-sp|--subsampling-prop)
			sub_prop=$2
			shift
			;;
		-mbs|--min-bases-sequenced)
			mbs=$2
			shift
			;;
		-mrc|--min-read-count)
			mrc=$2
			shift
			;;
		-p|--fb-ploidy)
			ploidy=$2
			shift
			;;
		-q|--fb-min-base-quality)
			minqual=$2
			shift
			;;
		-Q|--fb-mismatch-base-quality-threshold)
			mismatchqual=$2
			shift
			;;
		-m|--fb-min-mapping-quality)
			mapq=$2
			shift
			;;
		-mc|--fb-min-coverage)
			fb_min_cov=$2
			shift
			;;
		-n|--fb-best-alleles)
			nbest=$2
			shift
			;;
		-F|--fb-min-alternate-fraction)
			altfrac=$2
			shift
			;;
		-minlen|--min-locus-length)
			minlen=$2
			shift
			;;
		-miss|--max-missing)
			vcfmiss=$2
			shift
			;;
		*) echo "Unknown parameter passed: $1"
			exit 1
			;;
	esac
	shift
done

depends="makeblastdb blastn bedtools bwa samtools freebayes parallel vcf2fasta vcftools muscle Rscript AMAS.py awk" #samtools needs to be at least 1.10!

if [[ ${#ref} -le 1 && ${#popmap} -le 1 ]]; then
	echo "
	This help menu is shown because you did not specify the reference genome and population map, the two mandatory parameters for this pipeline.
	-r --reference-genome
							The file that contains the reference genome that should be used for both aligning the reads and haplotype calling.
	-popmap --population-map
							Tab-delimited text file with two columns: 1. identifier of samples without the file extension; 2. population identifier. The names in the first column must exactly match the file names that contain the reads. The population identifier will be used by freebayes. Please make sure IDs do not contain dots. This file can control which samples are being analysed.
	
	Other paramters that can be specified:
	-s --samples-directory
							The directory where fastq files are stored. The extension of the files can be .fq, .fastq,. fq.gz and fastq.gz. When paired-end reads are used read pairs must be indicated by .1 and .2; e.g. SAMPLE_ID.1.fq.gz and SAMPLE_ID.2.fq.gz, which is the default of Stacks' process_radtags. [default is the current directory]

	-o --output-directory
							Output directory to store all output files. Use absolute path. [default is the current directory]

	-mask --mask-reference 
							Valid options are "yes" or "no" (without quotes). Indicates if the inverted repeat of the chloroplast should be masked with N. Automatically removes only the largest duplication. [default no]

	-align --align-reads
							Valid options are "yes" or "no" (without quotes). Indicates if reads stored in --samples-directory and named by first column of --population-map should be aligned to reference. [default yes]						

	-np --number-of-processors
							The number of processors (threads) that will be used. [default 1]

	-type --sequencing-type						
							Valid options are SE (single-end) or PE (paired-end). In case of PE read pairs should be separated to SAMPLE_ID.1.fq.gz and SAMPLE_ID.2.fq.gz. [default SE]

	-bwa_k --min-seed-length
							Minimum seed length for the Burrows-Wheeler aligner. [default 19]

	-bwa_A --matching-score
							Matching-score for the Burrows-Wheeler aligner. [default 1]

	-bwa_B --mismatch-penalty
							Mismatch-penalty for the Burrows-Wheeler aligner. [default 4]

	-bwa_O --gap-open-penalty
							Gap-open penalty for the Burrows-Wheeler aligner. [default 6]

	-call --call-haplotypes
							Valid options are "yes" or "no" (without quotes). Indicates if haplotypes should be extracted to vcf and fasta format. The fasta format contains the entire sequences of the loci. It is suggested to first align reads to the reference, then investigate the read depth of loci. The bed files produced after alignment also report the bed coverage. 

	-mbc --min-bed-coverage 
							Minimal coverage of a locus to be included in the output haplotypes. If the bed coverage of a locus is higher in any sample than this value, it will be a subject of analysis in all the samples. As organellar DNA can be overrepresented, it can be a really high number (e.g. 500). When trying different thresholds it is not needed to align the reads to the reference again, only the genotype call step should be redone. [default 3]  

	-sp --subsampling-prop
							Float value between 0 and 1. To use less memory when calling haplotypes, each bed locus can be downsampled using samtools view by this proportion. 0.1 means: use 10% percent of all the reads found in a bed locus. The effect of downsampling is not tested properly, so use at your own risk and double check the results. [default 1.0]

	-mbs --min-bases-sequenced
							Defines the minimum number of bases sequenced in a sample to be included in the genotype calling. [default 0]

	-mrc --min-read-count
							Defines the minimum number of reads needed in a sample to be included in the genotype calling. [default 1]

	-p --fb-ploidy
							Ploidy level used by freebayes to call haplotypes. For organellar DNA it should be 1. [default 1]

	-q --fb-min-base-quality
							Minimal base quality to include in the analysis. [default 20]

	-Q --fb-mismatch-base-quality-threshold
							Minimal quality of a mismatched base to be accepted. [default 15]

	-m --fb-min-mapping-quality
							Minimal mapping quality for a read to be included in haplotype calling. [default 30]

	-F --fb-min-alternate-fraction
							Minimal supporting fraction of reads for an alternate allele. [default 0.4]

	-mc --fb-min-coverage
							Minimal coverage to accept when calling haplotypes. [default 5]

	-n --fb-best-alleles
							The number of most probable alleles to consider. [default 5]

	-minlen --min-locus-length
							Minimal length of loci to include in the final dataset. [default 100]

	-miss --max-missing
							Maximal missingness when filtering sites. [default 0.8]

	Dependencies are echo "$depends"
	Please make sure that samtools version is at least 1.10. It is not checked automatically.

	example:
	RADOrgMiner.sh -o ~/dataset/someoutput --mask-reference yes -r ~/dataset/reference_chloroplast.fa -align yes -call yes -np 24 -popmap ~/dataset/popmap -type PE -mbc 1000 -mbs 5000 -sp 0.1 -s ~/dataset/raw_reads
	"
	exit 1
fi

for i in $depends
do 
	if which $i; then
		echo $i found
	else 
		echo $i is not found
		echo "Please install $i or specify it in the '$PATH'"
		echo "The pipeline will continue now, but unexpected behaviour may follow"
	fi
done

if [[ ${#ref} -le 1 ]]; then
	echo "Please specify reference genome"
	exit 1
fi

if [[ ${#popmap} -le 1 ]]; then
	echo "Please specify popmap"
	exit 1
fi

#check if input files and refgenom exist

if [[ ! -d ${outdir} ]]; then
	echo "Output directory does not exist"
	exit 1
fi

echo "Output directory is $outdir" 
echo "Reference genome is $ref"
echo "Detection of inverted repeat is set to $mask"
echo "$np threads will be used"
echo "$popmap will be used for population map"
echo "Genotype calling is set to $call"
echo "Sequencig type is $type"
echo "BWA parameters are k=$bwa_k A=$bwa_A B=$bwa_B O=$bwa_O"
echo "Minimal bed coverage is $min_bed_cov"
echo "Subsampling proportion is $sub_prop"
echo "Ploidy is set to $ploidy"
echo "Minimal base quality for genotype calling is $minqual"
echo "Minimal base quality for genotype calling when a mismatch occurs is $mismatchqual"
echo "Minimal mapping quality for genotype calling is $mapq"
echo "Minimal coverage for genotype calling is $fb_min_cov"
echo "The best $nbest alleles will be considered"
echo "Minimal fraction of reads to support an alternate allele is $altfrac"
echo "Final missingness is set to $vcfmiss"
echo "Minimal number of basepairs required to be sequenced in a samples is set to $mbs"
echo "Minimal number reads required to be aligned in a samples is set to $mrc"
echo "Minimal bed coverage to consider a locus is $min_bed_cov"

echo "popmap is:"
cat $popmap

inds=`cut -f 1 $popmap | sort`

#depends: makeblastdb, blastn, bedtools
if [[ "$mask" == "yes" ]]; then													 #ALWAYS SPECIFY THE ORIGINAL, UNMASKED REFERENCE
	ref1=`echo $ref | cut -f 1 -d "."`
	if [ -f "${ref1}_nr.fasta" ]; then
	
		echo "Masked reference file found"
	
		echo "Assign ${ref1}_nr.fasta to reference database"
	
		ref_db=${ref1}_nr.fasta
	
	else
		echo "Attempting to find boundaries of inverted repeat in $ref then masking reference fasta file"
	
		makeblastdb -in $ref -dbtype nucl -out ${ref1}_blastdb
	
		blastn -db ${ref1}_blastdb -query $ref -out ${ref1}_selfblast.fmt6 -outfmt '6 qseqid qstart qend length evalue qseq pident' -num_threads $np
	
		cut -f 1-4 ${ref1}_selfblast.fmt6 | awk 'NR == 2{print;exit}' | cut -f 1-3 > ${ref1}_IR
	
		mv ${ref1}_IR ${ref1}_IR_boundary
	
		bedtools maskfasta -fi $ref -bed ${ref1}_IR_boundary -fo ${ref1}_nr.fasta 
	
		echo "Assigning ${ref1}_nr.fasta to reference database"
	
		ref_db=${ref1}_nr.fasta
	fi
elif [[ "$mask" == "no" ]]; then
	ref1=`echo $ref | cut -f 1 -d "."`
	if [ -f "${ref1}_nr.fasta" ]; then
	
		echo "Masked reference file found"
	
		echo "Assign ${ref1}_nr.fasta to reference database"
	
		ref_db=${ref1}_nr.fasta
	else
		echo "Masking of reference was not requested"
	
		echo "Assign $ref to reference database"
	
		ref_db=$ref
	fi
fi
#finds _nr.fasta if exists even if mask is not requested

if [[ "$map" == "yes" && ${#indir} -le 1 ]]; then
	echo "Please specify input directory of reads to be aligned"
	exit 1
fi

if [[ "$map" == "yes" && "$type" == "SE" ]]; then
	if [[ -d "${outdir}/map_reads" ]]; then
		rm -r ${outdir}/map_reads
	fi
	if [[ ! -d "${outdir}/unaligned" ]]; then
		mkdir ${outdir}/unaligned
	fi
	if [[ ! -d "${outdir}/aligned" ]]; then
		mkdir ${outdir}/aligned
	fi
	mkdir ${outdir}/map_reads
	
	bwa index $ref_db
	
	for i in $inds;
	do
		r1=`find $indir -maxdepth 1 -name ${i}.fq -or -name ${i}.fq.gz -or -name ${i}.fastq -or -name ${i}.fastq.gz`  #extension can be fq,fq.gz,fastq,fastq.gz
	
		if [[ ${#r1} -le 1 ]]; then
			echo "Sample $i was not found"
			exit 1
		fi
		
		echo "Aligning $r1 to $ref_db using bwa" 					      #sequence name must match with popmap

		bwa mem -t $np -k $bwa_k -A $bwa_A -B $bwa_B -O $bwa_O -R "@RG\tID:$i\tSM:$i\tPL:Illumina" $ref_db $r1 2> /dev/null |\
		samtools view -h -b -u -@ $np |\
		samtools sort -@ $np > ${outdir}/map_reads/${i}.bam

		samtools view -h -b -F 4 -@ $np ${outdir}/map_reads/${i}.bam > ${outdir}/aligned/${i}.bam
		samtools index -@ $np ${outdir}/aligned/${i}.bam
		#if bam is empty stop

		samtools view -f 4 -@ $np ${outdir}/map_reads/${i}.bam |\
		samtools sort -n -@ $np |\
		samtools fastq -@ $np -0 ${outdir}/unaligned/${i}.fq.gz 2> /dev/null

		samtools coverage ${outdir}/aligned/${i}.bam -m

		samtools depth -d 0 -a ${outdir}/aligned/${i}.bam > ${outdir}/aligned/${i}.depth

		echo "Per site read depth can be found in ${outdir}/aligned/${i}.depth"
	done

	for i in $inds
	do
		echo "
		bedtools bamtobed -i ${outdir}/aligned/${i}.bam |\
		cut -f 1-3 |\
		sort -n -k 2 |\
		uniq |\
		bedtools coverage -b ${outdir}/aligned/${i}.bam -a stdin -counts > ${outdir}/aligned/${i}.bed"
	done | parallel -j $np

	rm -r ${outdir}/map_reads

	echo "Depth of bed loci can be found in ${outdir}/aligned/SAMPLE_ID.bed"
	echo "Reads that did not align to the reference are stored in ${outdir}/unaligned"

	cd ${outdir}/aligned
	for i in *depth
	do
		echo $i > temp_${i}
		cut -f 3 $i >> temp_${i}
	done
	paste temp* > ind_depths.tsv
	rm temp*
	for i in $inds
	do
	 	echo ${i}.bam
	 	samtools coverage ${i}.bam
	done | grep -v ^# | perl -pe 's/.bam\n/\t/' > coverage_table.tsv


	echo "#!"$( which Rscript )"
	wd<-\"$( echo ${outdir}/aligned )\"
	setwd(wd)
	x<-read.table(\"ind_depths.tsv\", sep=\"\t\", header=T)
	col_names<-colnames(x)
	print(col_names)
	for (i in col_names){
	name<-paste(\"read_depth\", i, \".png\", sep=\"\")
	png(name, width=3000, height=1000, unit=\"px\")
	par(mar=rep(4,4))
	plot(x[,i], main=i, ylab=\"Read depth\", xlab=\"Position\", type=\"b\")
	dev.off()
	}
	means<-rowMeans(x)
	png(\"meandepth.png\", width=3000, height=1000, unit=\"px\")
	plot(means, type=\"b\", lty=1, lwd=2, cex=1, xlab=\"Position\", ylab = \"Mean read depth\")
	dev.off()
	" > plot_depth.R

	echo "Read depth of each site of samples can be checked by running "Rscript ${outdir}/aligned/plot_depth.R""
	echo "Coverage statistics can be found in ${outdir}/aligned/coverage_table.tsv"

	cd $cdir

elif [[ "$map" == "yes" && "$type" == "PE" ]]; then
	if [[ -d "${outdir}/map_reads" ]]; then
		rm -r ${outdir}/map_reads
	fi
	if [[ ! -d "${outdir}/unaligned" ]]; then
		mkdir ${outdir}/unaligned
	fi
	if [[ ! -d "${outdir}/aligned" ]]; then
		mkdir ${outdir}/aligned
	fi

	mkdir ${outdir}/map_reads

	bwa index $ref_db

	samtools faidx $ref_db

	for i in $inds;
	do
		r1=`find $indir -maxdepth 1 -name ${i}.1.fq -or -name ${i}.1.fq.gz -or -name ${i}.1.fastq -or -name ${i}.1.fastq.gz`
		r2=`find $indir -maxdepth 1 -name ${i}.2.fq -or -name ${i}.2.fq.gz -or -name ${i}.2.fastq -or -name ${i}.2.fastq.gz`
			
		if [[ ${#r1} -le 1 ]]; then
			echo "Sample $i was not found"
			exit 1
		fi

		echo "Aligning ${i}.1 & ${i}.2 to $ref_db" 	       #sequence name must match with popmap, but must end with 1.fq.gz and 2.fq.gz

		bwa mem -t $np -k $bwa_k -A $bwa_A -B $bwa_B -O $bwa_O -R "@RG\tID:$i\tSM:$i\tPL:Illumina" $ref_db $r1 $r2 2> /dev/null |\
		samtools view -h -b -u -@ $np |\
		samtools sort --threads $np > ${outdir}/map_reads/${i}.bam

		samtools view -h -b -F 12 -@ $np ${outdir}/map_reads/${i}.bam > ${outdir}/aligned/${i}.bam
		samtools index -@ $np ${outdir}/aligned/${i}.bam
		#if bam is empty stop

		samtools view -f 12 -@ np ${outdir}/map_reads/${i}.bam |\
		samtools sort -n -@ $np |\
		samtools fastq -@ $np -1 ${outdir}/unaligned/${i}.1.fq.gz -2 ${outdir}/unaligned/${i}.2.fq.gz 2> /dev/null
		#to emit singlets use samtools view -f 4 instead of -f 12
		nsites=`samtools depth ${outdir}/aligned/${i}.bam | awk '$3 > 0' | wc -l`
		echo "There are $nsites sites in the samples $i with read depth > 0"
		
		samtools coverage ${outdir}/aligned/${i}.bam -m

		samtools depth -d 0 -a ${outdir}/aligned/${i}.bam > ${outdir}/aligned/${i}.depth

		echo "Per site read depth can be found in ${outdir}/aligned/${i}.depth"
	done

	echo "Exporting read depth of bed loci"

	for i in $inds
	do
		echo "
		bedtools bamtobed -i ${outdir}/aligned/${i}.bam |\
		cut -f 1-3 |\
		sort -n -k 2 |\
		uniq |\
		bedtools coverage -b ${outdir}/aligned/${i}.bam -a stdin -counts > ${outdir}/aligned/${i}.bed"
	done | parallel -j $np

	rm -r ${outdir}/map_reads

	echo "Depth of bed loci can be found in ${outdir}/aligned/SAMPLE_ID.bed"
	echo "Reads that did not align to the reference are stored in ${outdir}/unaligned"

	cd ${outdir}/aligned
	for i in *depth
	do
		echo $i > temp_${i}
		cut -f 3 $i >> temp_${i}
	done
	paste temp* > ind_depths.tsv
	rm temp*
	for i in $inds
	do
	 	echo ${i}.bam
	 	samtools coverage ${i}.bam
	done | grep -v ^# | perl -pe 's/.bam\n/\t/' > coverage_table.tsv


	echo "#!"$( which Rscript )"
	wd<-\"$( echo ${outdir}/aligned )\"
	setwd(wd)
	x<-read.table(\"ind_depths.tsv\", sep=\"\t\", header=T)
	col_names<-colnames(x)
	print(col_names)
	for (i in col_names){
	name<-paste(\"read_depth\", i, \".png\", sep=\"\")
	png(name, width=3000, height=1000, unit=\"px\")
	par(mar=rep(4,4))
	plot(x[,i], main=i, ylab=\"Read depth\", xlab=\"Position\", type=\"b\")
	dev.off()
	}
	means<-rowMeans(x)
	png(\"meandepth.png\", width=3000, height=1000, unit=\"px\")
	plot(means, type=\"b\", lty=1, lwd=2, cex=1, xlab=\"Position\", ylab = \"Mean read depth\")
	dev.off()
	" > plot_depth.R

	echo "Read depth of each site of samples can be checked by running "Rscript ${outdir}/aligned/plot_depth.R""
	echo "Coverage statistics can be found in ${outdir}/aligned/coverage_table.tsv"

	cd $cdir

else

	echo "Aligning to reference was not requested. If you intended to map the reads please check if a typo happened when specifying --map-reads or --sequencing-type"

fi

if [[ $call == "yes" ]]; then

	grep -F -f <(awk -v x=$mbs -v y=$mrc '$6 > x && $5 > y' ${outdir}/aligned/coverage_table.tsv | cut -f 1) $popmap > ${outdir}/aligned/popmap_mbs

	popmap_mbs="${outdir}/aligned/popmap_mbs"
	inds_mbs=`cut -f 1 $popmap_mbs`
	echo "The popmap reduced by -mbs and -mrc is:"
	cat $popmap_mbs

	echo "Subsetting bed loci with a subsampling proportion of $sub_prop"
	for i in $inds_mbs
	do
		awk -v mincov="$min_bed_cov" '$4 > mincov' ${outdir}/aligned/${i}.bed |\
		bedtools merge |\
		bedtools coverage -b ${outdir}/aligned/${i}.bam -a stdin -counts > ${outdir}/aligned/${i}_merge.bedcov
	done

	cut -f 1-3 ${outdir}/aligned/*_merge.bedcov | bedtools sort | bedtools merge > ${outdir}/aligned/merged_alignments.bed
	
	echo "bed regions with a minimum bed coverage of $min_bed_cov are:"
	cat ${outdir}/aligned/merged_alignments.bed

	bedtools getfasta -fi $ref_db -fo ${outdir}/aligned/reference_subset.fa -bed ${outdir}/aligned/merged_alignments.bed
	bwa index ${outdir}/aligned/reference_subset.fa

	proplarger=`awk -v x=$sub_prop 'BEGIN { print (x > 1.0) ? "yes" : "no" }'`

	if [[ $proplarger == "yes" ]]; then
		echo "Subsampling of reads can not be done with a proportion larger than 1.0"
		exit 1
	fi

	for i in $inds_mbs
	do
		echo $i
		cat ${outdir}/aligned/merged_alignments.bed |\
		while read line
		do 
		coord=`echo $line | cut -f 1-3 -d " "`
		#echo $coord
		samtools view -h -b -L <(echo $coord) -s $sub_prop ${outdir}/aligned/${i}.bam > "${outdir}/aligned/ssmp_${i}_${sub_prop}_${coord}.bam"
		done
			
		samtools merge -b <(ls ${outdir}/aligned/ssmp_${i}*.bam) -f ${outdir}/aligned/${i}_subsampled.bam
		
		rm ${outdir}/aligned/ssmp_${i}_*.bam
			
		samtools index -@ $np ${outdir}/aligned/${i}_subsampled.bam
	done

	for i in $inds_mbs
	do
		bwa mem -t $np -R "@RG\tID:$i\tSM:$i\tPL:Illumina" ${outdir}/aligned/reference_subset.fa <(samtools bam2fq -@ $np ${outdir}/aligned/${i}_subsampled.bam) 2> /dev/null | samtools view -h -b -@ $np | samtools sort -@ $np > ${outdir}/aligned/${i}_subset.bam
		samtools index -@ $np ${outdir}/aligned/${i}_subset.bam
		bedtools bamtobed -i ${outdir}/aligned/${i}_subset.bam | bedtools merge > ${outdir}/aligned/${i}_subset.bed
	done

	samtools faidx ${outdir}/aligned/reference_subset.fa
	awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${outdir}/aligned/reference_subset.fa.fai > ${outdir}/aligned/subset_alignments.bed

	rm ${outdir}/aligned/*subsampled.bam
	rm ${outdir}/aligned/*subsampled.bam.bai

	ls ${outdir}/aligned/*_subset.bam | grep "$inds_mbs" > ${outdir}/aligned/bamlist

	if [[ ! -d ${outdir}/aligned/bed_loci ]]; then
		mkdir ${outdir}/aligned/bed_loci
	elif [[ -d ${outdir}/aligned/bed_loci ]]; then
		rm ${outdir}/aligned/bed_loci/*
	fi

	cat ${outdir}/aligned/subset_alignments.bed |\
	while read line
	do 
		echo "$line" > "${outdir}/aligned/bed_loci/loc_${line}"
	done

	cd ${outdir}/aligned/bed_loci
	rename 's/\t/_/' ./* 
	rename 's/\t/-/' ./* 
	rename 's/..:/_/' ./*
	loc=`ls`
	cd $cdir

	if [[ ! -d ${outdir}/aligned/vcf_loci ]]; then
		mkdir ${outdir}/aligned/vcf_loci
	elif [[ -d ${outdir}/aligned/vcf_loci ]]; then
		rm ${outdir}/aligned/vcf_loci/*
	fi

	echo "Calling haplotypes with freebayes"

	for i in $loc
	do
		echo "freebayes -f ${outdir}/aligned/reference_subset.fa -L ${outdir}/aligned/bamlist --ploidy $ploidy -q $minqual -Q $mismatchqual -m $mapq --min-coverage $fb_min_cov -w -j -V -E -1 -n $nbest -F $altfrac --populations $popmap_mbs --report-monomorphic -t ${outdir}/aligned/bed_loci/${i} > ${outdir}/aligned/vcf_loci/${i}.vcf"
	done | parallel -j $np

	if [[ ! -d ${outdir}/fasta_loci ]]; then
		mkdir ${outdir}/fasta_loci
	elif [[ -d ${outdir}/fasta_loci ]]; then
		rm ${outdir}/fasta_loci/* #by using find it would not run into an error when dir. is empty
	fi

	ref_id=`head -n 1 $ref_db | sed 's/>//' | cut -f 1 -d " "`
	cd ${outdir}/aligned/vcf_loci

	echo "Exporting sequnces of bed loci to fasta and variant sites to vcf"
	
	ncomment=`for i in *.vcf; do grep "^#" $i | wc -l; done | sort | head -n 1`

	find ./ -name "*.vcf" -type f -exec awk -v x=$ncomment 'NR==x+1{exit 1}' {} \; -exec echo rm {} \; > rmvcf

	if [[ $(wc -l <rmvcf) -ge 1 ]]; then 
		sh rmvcf &> /dev/null
		rm rmvcf
	else
		rm rmvcf
	fi

	vcfloc=`ls *.vcf | sed 's/.vcf//'`

	for i in $vcfloc
	do
		#echo $i
		vcftools --vcf ${i}.vcf --max-missing $vcfmiss --recode --out $i 2> /dev/null
	done

	nrecode_comment=`for i in *recode.vcf
	do
		grep "^#" $i | wc -l
	done | sort | head -n 1`

	find ./ -name "*.recode.vcf" -type f -exec awk -v x=$nrecode_comment -v y=$minlen 'NR==x+y{exit 1}' {} \; -exec echo rm {} \; > rmrecodevcf
	
	if [[ $(wc -l <rmrecodevcf) -ge 1 ]]; then 
		sh rmrecodevcf &> /dev/null
		rm rmrecodevcf
	else	
		rm rmrecodevcf
	fi

	recodeloc=`ls *.recode.vcf | sed 's/.recode.vcf//'`
	for i in $recodeloc
	do
		vcf2fasta -f ${outdir}/aligned/reference_subset.fa -p $i -P 1 -n N ${outdir}/aligned/vcf_loci/${i}.recode.vcf
		cat ${i}*.fa | perl -pe "s/_"${ref_id}".*//" | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > ${outdir}/fasta_loci/unaligned_${i}.fa

		for l in $inds_mbs
		do
		grep -v "^#" ${outdir}/aligned/vcf_loci/${i}.recode.vcf | cut -f 2 | sed -n '1p;$p' | perl -pe 's/\n/\t/' 
		echo 
		done | sed 's/1\t/0\t/' > ${outdir}/aligned/vcf_loci/${i}.range

		paste <(echo $inds_mbs | perl -pe 's/ /\n/g') <(cat ${i}.range) > ${outdir}/fasta_loci/${i}.bedlocus
	done

	rm *range
	rm *fa

	cd ${outdir}/fasta_loci

	empty_loci=`wc -l *fa | awk '$1 == 0' | sed 's/^ *//' | cut -f 2 -d" "`
	if [[ ${#empty_loci} -ge 1 ]]; then
		echo $empty_loci | xargs rm
	fi

	for i in $recodeloc
	do 
		muscle -in unaligned_${i}.fa -out break_${i}.fa &> /dev/null
		awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' break_${i}.fa > miss_${i}.fa
		bedtools getfasta -fi miss_${i}.fa -fo ${i}.fa -bed ${i}.bedlocus &> /dev/null
		sed -i 's/:.*//' ${i}.fa
		sed -i -e '/^[^>]/s/N/-/g' ${i}.fa #missing is coded as "-"
	done
	rm unaligned_*.fa
	rm break_*.fa
	rm miss_*.fa
	rm *.fai
	rm *bedlocus

	AMAS.py summary -i loc*fa -f fasta -d dna

	cd $outdir
	AMAS.py concat -i fasta_loci/loc*fa -f fasta -d dna -t concat_loci.fa -p concat_loci.parts
	sed -i 's/^/DNA,/' concat_loci.parts

	cd $cdir

	echo "Sequence of each locus can be found in ${outdir}/fasta_loci"
	echo "Concatenated sequence of all loci is written to ${outdir}/concat_loci.fa"
	echo "Summary of alignments can be found in ${outdir}/fasta_loci/summary.txt"
fi


echo "Run ended at $(date)"

#end
