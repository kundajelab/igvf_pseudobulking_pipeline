#!/bin/bash

basedir="${1}"

# get only final fragments (overwrite nonsorted by sorted)
echo -e "\t\t- renaming fragment files"
for f in ${basedir}/pseudobulked_fragments/*-sorted.tsv; do
	name=$(basename "$f" -sorted.tsv)
	mv "$f" "${basedir}/pseudobulked_fragments/${name}.tsv"
done

# get the final peakset+signal track
echo -e "\t\t- renaming peak files"
mkdir "${basedir}/peaks_temp"
mkdir "${basedir}/peaks_bw_temp"
for f in ${basedir}/peaks/*-peaks_overlap_filtered.narrowPeak; do
	name=$(basename "$f" -peaks_overlap_filtered.narrowPeak)
	mv "$f" "${basedir}/peaks_temp/${name}.narrowPeak"
done
for f in ${basedir}/peaks/*-pval.bw; do
	mv "$f" "${basedir}/peaks_bw_temp/"
done
rm ${basedir}/peaks/*
for f in ${basedir}/peaks_temp/*; do
	mv "$f" "${basedir}/peaks/"
done
for f in ${basedir}/peaks_bw_temp/*; do
	mv "$f" "${basedir}/peaks/"
done
rm -r "${basedir}/peaks_temp"
rm -r "${basedir}/peaks_bw_temp"

