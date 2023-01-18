#!/bin/bash
get_input () {
        # Function for doing getopts

        # Setting default number of threads as 4
        threads=4

        while getopts "t:q:p:u:o:m:vh" option
        do
                case $option in
                        t) threads=$OPTARG;;
			q) quastPath=$OPTARG;;
			p) pathToInputFiles=$OPTARG;;
                        u) pathToContig=$OPTARG;;
                        o) outputFolder=$OPTARG;;
			h) info_usage=1;;
                        v) verbose=1;;
                        *) echo "Incorrect arguments used. Use -h to know about the arguments available."
                esac
        done
	
        if ((info_usage)); then
                echo -e "The script contains a pipeline for running quast evaluation.\nRun the script in the following format after giving the script executable permission:\n./run_quast.sh\nArguments available\n\t-q <Path to quast.py file>\n\t-p <Path to input files>\n\t-u <Path to the contig files>\n\t-o <Path to output folder>\n\t-t <Number of threads>\n\t-v \tVerbose mode\n\t-h\tPrint usage information"
                exit
        fi
}

get_input "$@"

if ((verbose)); then
	echo "\nMaking output directory\n"
fi

if ((verbose)); then
        echo "\nStarted process of identification of Assembler to use\n"
fi

varn=""
varf=""

mkdir ${outputFolder}/quast_output

for name in `ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw` 
do 
	varn+=" ${name}"
	varf+=" ${outputFolder}/assembled_outputs/${name}_assembled.fasta"
done

varn=`echo $varn | tr ' ' ', '`

${quastPath} -o ${outputFolder}/quast_output -t ${threads} -l ${varn} ${varf}

mv ${outputFolder}/quast_output/report.txt ${outputFolder}/assembly_statistics.txt
rm -r ${outputFolder}/quast_output
