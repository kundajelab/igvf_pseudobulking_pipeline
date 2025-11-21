#!/bin/bash

basedir="${1}"

mkdir -p "${basedir}/separated_pseudorep1"
mkdir -p "${basedir}/separated_pseudorep2"
mkdir -p "${basedir}/separated_pseudorepT"
mkdir -p "${basedir}/separated_fragments"

mkdir -p "${basedir}/pseudobulked_pseudorep1"
mkdir -p "${basedir}/pseudobulked_pseudorep2"
mkdir -p "${basedir}/pseudobulked_pseudorepT"
mkdir -p "${basedir}/pseudobulked_fragments"

mkdir -p "${basedir}/peaks"

mkdir -p "${basedir}/pseudobulked_rna"

mkdir -p "${basedir}/rna_qc_reports"
mkdir -p "${basedir}/atac_qc_reports"
