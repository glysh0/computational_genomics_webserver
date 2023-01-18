#!/usr/bin/bash
#activate environment
#source activate geneprediction
v=0
g=0
r=0
t=0
while getopts ":hi:o:grtv" opt
do
    case $opt in
       h )
         echo "run_dfast.sh -i <input of file> -o <path to output> -g <run GMS2 instead of prodigal> -r <run RNAmmer instead of barrnarp> -t <run tRNAscan-SE instead of Aragorn> -v <verbose> "
         exit;;
       i )
         i=$OPTARG
         ;;
       o )
         o=$OPTARG
         ;;
       r )
         r=1
         ;;
       g )
         g=1
         ;;
       t )
         t=1
         ;;
       v )
         v=1
         ;;
     esac
done
genome=$i
output=$o
#get file name
genome_name=$( echo $genome | sed -e "s|.*\/||g" -e "s|.fasta||g" )

# create output directory if doesn't exist
if [ ! -d ${output}amino_acids ]
then
    echo $output
    mkdir -p ${output}amino_acids
fi

if [ "$v" -eq 1 ]
then
    echo Analyzing $genome_name
fi

# set CDS prediction tool
if [ "$g" -eq 1 ]
then
    cds_tool='--use_genemarks2 bact'
else
    cds_tool='--use_prodigal'
fi

base_command="/projects/VirtualHost/predicta/html/miniconda3/bin/python3.7 /projects/VirtualHost/predicta/html/Team1-WebServer/tools/dfast_core/dfast -g ${genome} -o tmp_dfast --organism Ecoli --minimum_length 120 ${cds_tool} --gcode 11 --cpu 8"

# run dfast with RNAmmer instead of barrnap
if [ "$r" -eq 1 ]
then
    base_command="${base_command} --use_rnammer bact"
    echo $base_command
fi

# run dfast with trnascane instead of aragorn
if [ "$t" -eq 1 ]
then
    base_command="${base_command} --use_trnascan bact"
fi

$base_command

# move files in final path and clean up
mv tmp_dfast/genome.gff ${output}${genome_name}.gff
mv tmp_dfast/protein.faa ${output}amino_acids/${genome_name}_protein.faa
mv tmp_dfast/cds.fna ${output}${genome_name}_cds.fna
mv tmp_dfast/rna.fna ${output}${genome_name}_rna.fna
rm -r tmp_dfast
#rm -r gene_predicted_files
