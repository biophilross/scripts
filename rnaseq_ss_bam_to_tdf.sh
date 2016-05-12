#!/bin/bash

set -ueo pipefail

if [ $# -ne 6 ]; then
	echo -e "\nusage: $(basename $0) <in.bam> <genome.windows> <minimum.mapq> <genome.lengths> <genome.fasta> <outdir>\n"
	echo -e "This script will take a paired-end strand specific RNA-Seq BAM file"
	echo -e "and convert it to TDF format to visualize in IGV by strand. You can choose a window size by"
	echo -e "using the bedtools makewindows command and using that as an input to <genome.windows>.\n"
	exit
fi

INBAM=$1
WINDOWS=$2
QUALITY=$3
GENOME_LENGTHS=$4
GENOME_FASTA=$5
OUTDIR=$6

# create output directory if it doesn't exist
if [ ! -d "${OUTDIR}" ]; then
	mkdir -p "${OUTDIR}"
fi

prefix=$(basename ${INBAM%.bam})

# count number of reads for normalizing coverage
num_reads=$(samtools view -q "${QUALITY}" -c "${INBAM}")

# split file into positive and negative strand
>&2 echo -e "Processing ${INBAM}...\n"
samtools view -q "${QUALITY}" -h "${INBAM}" | awk 'BEGIN {OFS=FS="\t"} /^@/ || /XS:A:\+/{print $0}' | samtools view -b - | \
	bedtools coverage -a "${WINDOWS}" -b - -g "${GENOME_LENGTHS}" -mean -split -bed > "${OUTDIR}/${prefix}.pos.cov"
samtools view -q "${QUALITY}" -h "${INBAM}" | awk 'BEGIN {OFS=FS="\t"} /^@/ || /XS:A:\-/{print $0}' | samtools view -b - | \
	bedtools coverage -a "${WINDOWS}" -b - -g "${GENOME_LENGTHS}" -mean -split -bed > "${OUTDIR}/${prefix}.neg.cov"
>&2 echo -e "Finished separating ${INBAM} by strand and calculating coverage\n"

# convert to igv friendly format
>&2 echo -e "Converting ${INBAM} to igv friendly format...\n"
cat "${OUTDIR}/${prefix}.pos.cov" | awk -v nr="${num_reads}" 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,"cov",$4 *(1000000/nr)}' > "${OUTDIR}/${prefix}.pos.igv"
cat "${OUTDIR}/${prefix}.neg.cov" | awk -v nr="${num_reads}" 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,"cov",-1 * $4 * (1000000/nr)}' > "${OUTDIR}/${prefix}.neg.igv"
>&2 echo -e "Finished converting to igv friendly format\n"

# convert to tdf
>&2 echo -e "Converting ${INBAM} to TDF format using igvtools...\n"
igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/${prefix}.pos.igv" "${OUTDIR}/${prefix}_pos.tdf" "${GENOME_FASTA}"
igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/${prefix}.neg.igv" "${OUTDIR}/${prefix}_neg.tdf" "${GENOME_FASTA}"
>&2 echo -e "Finished converting ${INBAM} to TDF format\n"

# remove temporary files
>&2 echo -e  "Cleaning up..."
rm -rf "${OUTDIR}/${prefix}.pos.cov"
rm -rf "${OUTDIR}/${prefix}.neg.cov"
rm -rf "${OUTDIR}/${prefix}.pos.igv"
rm -rf "${OUTDIR}/${prefix}.neg.igv"
rm -rf "./igv.log"
>&2 echo -e "All done!"

exit
