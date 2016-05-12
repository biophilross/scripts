#!/bin/bash

set -ueo pipefail

if [ $# -ne 4 ]; then
	echo -e "\nusage: $(basename $0) <in.bam> <genome.lengths> <genome.fasta> <outdir>\n"
	echo -e "This script will take any paired-end alignment file"
	echo -e "and convert it to TDF format for visualize in IGV by strand\n"
	exit
fi

INBAM=$1
GENOME_LENGTHS=$2
GENOME_FASTA=$3
OUTDIR=$4

# create output directory if it doesn't exist
if [ ! -d "${OUTDIR}" ]; then
	mkdir -p "${OUTDIR}"
fi

prefix=$(basename ${INBAM%.bam})

# split file into positive and negative strand
>&2 echo -e "Processing ${INBAM}...\n"
samtools view -h -F 16 -b "${INBAM}" | \
	bedtools genomecov -ibam - -g "${GENOME_LENGTHS}" -d > "${OUTDIR}/${prefix}.pos.cov"
samtools view -h -f 16 -b "${INBAM}" | \
	bedtools genomecov -ibam - -g "${GENOME_LENGTHS}" -d > "${OUTDIR}/${prefix}.neg.cov"
>&2 echo -e "Finished separating ${INBAM} by strand and calculating coverage\n"

# convert to igv friendly format
>&2 echo -e "Converting ${INBAM} to igv friendly format...\n"
cat "${OUTDIR}/${prefix}.pos.cov" | awk 'BEGIN {OFS=FS="\t"} {print $1,$2-1,$2,"cov",$3}' > "${OUTDIR}/${prefix}.pos.igv"
cat "${OUTDIR}/${prefix}.neg.cov" | awk 'BEGIN {OFS=FS="\t"} {print $1,$2-1,$2,"cov","-"$3}' > "${OUTDIR}/${prefix}.neg.igv"
>&2 echo -e "Finished converting to igv friendly format\n"

# convert to tdf
>&2 echo -e "Converting ${INBAM} to TDF format using igvtools...\n"
igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/${prefix}.pos.igv" "${OUTDIR}/${prefix}.pos.tdf" "${GENOME_FASTA}"
igvtools toTDF -z 5 -f min,max,mean "${OUTDIR}/${prefix}.neg.igv" "${OUTDIR}/${prefix}.neg.tdf" "${GENOME_FASTA}"
>&2 echo -e "Finished converting ${INBAM} to TDF format\n"

# remove temporary files
>&2 echo -e  "Cleaning up..."
rm -rf "${OUTDIR}/${prefix}.pos.cov"
rm -rf "${OUTDIR}/${prefix}.neg.cov"
rm -rf "${OUTDIR}/${prefix}.pos.igv"
rm -rf "${OUTDIR}/${prefix}.neg.igv"
rm -rf "./igv.log"
>&2 echo -e "All done!\n"

exit
