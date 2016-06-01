#!/usr/bin/env bash

set -euo pipefail

if [ $# -lt 3 ]; then
  echo -e "\nusage: $(basename $0) <ref.fasta> <read_length> <max_mismatches> [threads]\n"
  echo -e "Use this script to create a mappability track in bigWig format. Remember to" 
  echo -e "edit the chromosome names of your fasta file if need be to remove all instances"
  echo -e "of the '|' character. wigToBigWig doesn't like it..."
  exit
fi

FASTA=$1
READ_LENGTH=$2
MISMATCHES=$3
THREADS=${4:-1}

prefix="./$(basename ${FASTA%%.*})_k${READ_LENGTH}_m${MISMATCHES}"

echo -e "\nCreate appropriate FASTA file..."
tmpfile=$(mktemp ./$(basename $0).XXXXXX)
awk '{ if($0 ~ /^>/) {split($0, a, " | "); print a[1]} else {print $0} }' "${FASTA}" > "$tmpfile"

echo -e "\nIndexing FASTA file..."
gem_wrapperh.sh gem-indexer -i "${tmpfile}" -o "${prefix}"

echo -e "\nCalculating mappability..."
gem_wrapper.sh gem-mappability -I "${prefix}.gem" -o "${prefix}" -l "${READ_LENGTH}" -T "${THREADS}"

echo -e "\nConverting to wig format..."
gem_wrapper.sh gem-2-wig -I "${prefix}.gem" -i "${prefix}.mappability" -o "${prefix}"

echo -e "\nConverting to SAM format... "
gem_wrapper.sh gem-2-sam -I "${prefix}.gem" -i "${prefix}.mappability" -o "${prefix}"

echo -e "\nConverting to bigWig format..."
wigToBigWig "${prefix}.wig" "${prefix}.sizes" "${prefix}.bw"

echo -e "\nCleaning up..."
rm -rf "$tmpfile"
echo -e "\nDone!"
