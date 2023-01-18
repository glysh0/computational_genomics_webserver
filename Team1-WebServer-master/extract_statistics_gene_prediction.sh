#!/bin/bash
genome=$1
file_name=$( echo $genome | sed -e 's/\/.*\///g' )
# CDS
nr_cds=$( cat ${genome} | grep 'CDS' | wc -l )

# tRNA
nr_trna=$( cat ${genome} | grep 'tRNA' | wc -l )

# rRNA
nr_rrna=$( cat ${genome} | grep 'rRNA' | wc -l )

# CRISPR
nr_crispr=$( cat ${genome} | grep 'repeat_region' | wc -l )

echo ${file_name} ${nr_cds} ${nr_trna} ${nr_rrna} ${nr_crispr}
