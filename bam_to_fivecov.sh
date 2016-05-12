#!/bin/bash

set -ueo pipefail

if [ $# -ne 4 ]; then
	echo -e "\nusage: $(basename $0) <in.bam> <ext_length> <genome.lengths> <outdir>\n"
	echo -e "This script will take a single-end or paired-end ChIP-Seq BAM file"
	echo -e "separate it by strand, and create a bedgraph file of the 5' positions of"
	echo -e "reads extended out by <ext_length> for each strand\n"
	exit
fi

INBAM=$1
EXT_LENGTH=$2
GENOME_LENGTHS=$3
OUTDIR=$4

# create output directory if it doesn't exist
if [ ! -d "${OUTDIR}" ]; then
	mkdir -p "${OUTDIR}"
fi

prefix=$(basename ${INBAM%.*})

# split file into positive and negative strand
echo -e "Processing ${INBAM}...\n"
samtools view -F 16 "${INBAM}" | \
	awk 'BEGIN {OFS=FS="\t"} {print $3, $4, $4 + 1}' | \
	bedtools flank -i - -g "${GENOME_LENGTHS}" -l 0 -r "${EXT_LENGTH}" | \
	bedtools genomecov -i - -g "${GENOME_LENGTHS}" -bga > "${OUTDIR}/${prefix}_positive_5cov_ext${EXT_LENGTH}bps.bedgraph" 

samtools view -f 16 "${INBAM}" | \
	awk 'BEGIN {OFS=FS="\t"} {print $3, $4, $4 + 1}' | \
	bedtools flank -i - -g "${GENOME_LENGTHS}" -l 0 -r "${EXT_LENGTH}" | \
	bedtools genomecov -i - -g "${GENOME_LENGTHS}" -bga | \
	awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3, -$4}' > "${OUTDIR}/${prefix}_negative_5cov_ext${EXT_LENGTH}bps.bedgraph" 
echo -e "Finished separating ${INBAM} by strand and calculating coverage to creage bedgraph files\n"

# This code below I'm keeping there to upgrade to creadting TDF files for
# visualization in IGV. In the future I'll give an option as to whether I want to
# output bedgraph or TDF files

# convert to igv friendly format
#printf "Converting ${INBAM} to igv friendly format...\n\n"
#cat "${OUTDIR}/$prefix.pos.cov" | awk 'BEGIN {OFS=FS="\t"} {print $1,$2-1,$2,"cov",$3}' > "${OUTDIR}/$prefix.pos.igv"
#cat "${OUTDIR}/$prefix.neg.cov" | awk 'BEGIN {OFS=FS="\t"} {print $1,$2-1,$2,"cov",$3}' > "${OUTDIR}/$prefix.neg.igv"
#printf "Finished converting to igv friendly format\n\n"

# convert to tdf
#printf "Converting ${INBAM} to TDF format using igvtools...\n\n"
#igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/$prefix.pos.igv" "${OUTDIR}/$prefix.5.pos.tdf" "${GENOME_FASTA}"
#igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/$prefix.neg.igv" "${OUTDIR}/$prefix.5.neg.tdf" "${GENOME_FASTA}"
#printf "Finished converting ${INBAM} to TDF format\n\n"

# remove temporary files
#printf "Cleaning up...\n"
#rm -rf "${OUTDIR}/*.cov"
#rm -rf "${OUTDIR}/*.igv"
#rm -rf "igv.log"
echo -e "All done!\n"

exit
