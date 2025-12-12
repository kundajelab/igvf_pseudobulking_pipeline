#!/bin/bash

basedir="${1}"

rm -r "${basedir}/separated_pseudorep1"
rm -r "${basedir}/separated_pseudorep2"
rm -r "${basedir}/separated_pseudorepT"
rm -r "${basedir}/separated_fragments"
rm -r "${basedir}/atac_qc_reports"

rm -r "${basedir}/pseudobulked_pseudorep1"
rm -r "${basedir}/pseudobulked_pseudorep2"
rm -r "${basedir}/pseudobulked_pseudorepT"
rm -r "${basedir}/pseudobulked_fragments"

rm -r "${basedir}/peaks"

rm -r "${basedir}/pseudobulked_rna"
rm -r "${basedir}/rna_qc_reports"

rm "${basedir}/step5_complete.txt"