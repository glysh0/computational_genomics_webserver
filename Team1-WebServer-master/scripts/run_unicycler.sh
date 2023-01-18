#!/bin/bash

## Run genome assembly pipeline

get_input () {
        # Function for doing getopts

        # Setting default number of threads as 4
        threads=4

        while getopts "t:g:p:o:q:k:v:h" option
        do
                case $option in
                        t) threads=$OPTARG;;
                        p) pathToInputFiles=$OPTARG;;
                        g) genomeAssembler=$OPTARG;;
                        o) outputFolder=$OPTARG;;
                        k) kmersize=$OPTARG;;
			q) qualityControl=$OPTARG;;
                        h) info_usage=1;;
                        v) verbose=$OPTARG;;
                        *) echo "Incorrect arguments used. Use -h to know about the arguments available."
                esac
        done
        if ((info_usage)); then
                echo -e "The script contains a Genome Assembly pipeline.\nRun the script in the following format after giving the script executable permission:\n./run_unicycler.sh\nArguments available\n\t-p <Path to folder containing input fastq forward and backward reads annotated as <name>_1 and <name>_2 and are gzipped>\n\t-q\tSpecify 1 to perform quality control and trimming\n\t-k <Number of k-mer steps to use in assembly (default: 10)>\n\t-g <Specify genome assembler to use, with options as follows:\n\t\t1) a\tauto\n\t\t2) u\tUnicycler\n\t\t3) s\tSPAdes\n\t-o <Path to output folder>\n\t-t <Number of threads>\n\t-v \tVerbose mode\n\t-h\tPrint usage information"
                exit
        fi
}

check_files () {
        # Function for checking for presence of input files and output folder

        if ((verbose));then
                echo "Checking if input arguments are correct"
        fi
        if test -d "$pathToInputFiles"; then
                if ((verbose)); then
                        echo "Path: $pathToInputFiles exists"
                fi
                if [ "$(ls -A $pathToInputFiles | grep fq.gz)" ]; then
                        :
                else
                        echo "No fastq.gz files in $pathToInputFiles exist"
                        exit 1
                fi
        else
                echo "Path: $pathToInputFiles doesn't exist"
                exit 1
        fi
	        if ((verbose));then
                echo "Checking if input arguments are correct"
        fi
        if test -d "$outputFolder"; then
                if ((verbose)); then
                        echo "Path: $outputFolder exists"
                fi
        else
                echo "Path: $outputFolder doesn't exist"
                exit 1
        fi
	spadesPath=/projects/VirtualHost/predicta/html/miniconda3/bin/spades.py
	unicyclerPath=/projects/VirtualHost/predicta/html/miniconda3/bin/unicycler
	fastp_path=/projects/VirtualHost/predicta/html/miniconda3/bin/fastp
        if test -f "$unicyclerPath"; then
                if ((verbose)); then
                        echo "Unicycler path  exists"
                fi
        else
                echo "Unicycler path doesn't exist"
                exit 1
        fi
        if test -f "$spadesPath"; then
                if ((verbose)); then
                        echo "SPAdes path exists"
                fi
        else
                echo "SPAdes path doesn't exist"
                exit 1
        fi
        if test -f "$fastp_path"; then
                if ((verbose)); then
                        echo "FASTp path exists"
                fi
        else
                echo "FASTp path doesn't exist"
                exit 1
        fi
		
}



quality_control () {
        if ((verbose)); then
                echo "\nExecuting quality control and trimming\nCreating folders for output:"
        fi
        mkdir ${outputFolder}fastp_outputs
        mkdir ${outputFolder}trimmed_reads
        ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa mkdir ${outputFolder}fastp_outputs/gwa_report
        ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${fastp_path} -w ${threads} -i ${pathToInputFiles}/gwa_1.fq.gz -I ${pathToInputFiles}/gwa_2.fq.gz -5 -3 -c -o ${outputFolder}trimmed_reads/gwa_1.fq.gz --out2 ${outputFolder}trimmed_reads/gwa_2.fq.gz -j ${outputFolder}fastp_outputs/gwa_report/fastp.json -h ${outputFolder}fastp_outputs/gwa_report/fastp.html
}

genome_assembly () {
        if ((verbose)); then
                echo "\nExecuting genome assembly with tool "
        fi
	mkdir ${outputFolder}assembled_outputs 
        if [ "$genomeAssembler" == "s" ]; then
                echo "SPAdes"
		if ((qualityControl)); then
                        if ((verbose)); then
                                echo "\nStarted genome assembly\n"
                        fi
                        ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa mkdir ${outputFolder}gwa_output
                        if ((kmersize)); then
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${spadesPath} -1 ${outputFolder}trimmed_reads/gwa_1.fq.gz -2 ${outputFolder}trimmed_reads/gwa_2.fq.gz -k $kmersize -t ${threads} -o ${outputFolder}gwa_output
                        else
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${spadesPath} -1 ${outputFolder}trimmed_reads/gwa_1.fq.gz -2 ${outputFolder}trimmed_reads/gwa_2.fq.gz -t ${threads} -o ${outputFolder}gwa_output
			fi
			rm -r ${pathToInputFiles}trimmed_reads
                        for v in `ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw`
                        do
                                cp ${outputFolder}${v}_output/contigs.fasta ${outputFolder}assembled_outputs/${v}_assembled.fasta
				cp ${outputFolder}fastp_outputs/${v}_report/fastp.html ${outputFolder}${v}_fastp_report.html 
                                rm -r ${outputFolder}${v}_output
                        done
			rm -r ${outputFolder}fastp_outputs/
                else
                        if ((verbose)); then
                                echo "\nStarted genome assembly\n"
                        fi

                        ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa mkdir ${outputFolder}gwa_output
			if ((kmersize)); then
                        	ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${spadesPath} -1 ${pathToInputFiles}/gwa_1.fq.gz -2 ${pathToInputFiles}/gwa_2.fq.gz -k $kmersize -t ${threads} -o ${outputFolder}gwa_output
			else
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${spadesPath} -1 ${pathToInputFiles}/gwa_1.fq.gz -2 ${pathToInputFiles}/gwa_2.fq.gz -t ${threads} -o ${outputFolder}gwa_output
			fi
			rm -r ${pathToInputFiles}trimmed_reads
                        for v in `ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw`
                        do
                                cp ${outputFolder}${v}_output/contigs.fasta ${outputFolder}assembled_outputs/${v}_assembled.fasta
				rm -r ${outputFolder}${v}_output
                        done
                fi

        elif [ "$genomeAssembler" == "a" ]; then
                echo "Unicycler\n"
                if ((qualityControl)); then
                        if ((verbose)); then
                                echo "\nStarted genome assembly\n"
                        fi
			ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa mkdir ${outputFolder}gwa_output
			if ((kmersize)); then
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${unicyclerPath} -1 ${outputFolder}trimmed_reads/gwa_1.fq.gz -2 ${outputFolder}trimmed_reads/gwa_2.fq.gz --kmer_count $kmersize -t ${threads} -o ${outputFolder}gwa_output
			else
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${unicyclerPath} -1 ${outputFolder}trimmed_reads/gwa_1.fq.gz -2 ${outputFolder}trimmed_reads/gwa_2.fq.gz -t ${threads} -o ${outputFolder}gwa_output
			fi
			rm -r ${pathToInputFiles}trimmed_reads
			for v in `ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw`
			do
        			cp ${outputFolder}${v}_output/assembly.fasta ${outputFolder}assembled_outputs/${v}_assembled.fasta
        			cp ${outputFolder}fastp_outputs/${v}_report/fastp.html ${outputFolder}${v}_fastp_report.html
				rm -r ${outputFolder}${v}_output
			done
			rm -r ${outputFolder}fastp_outputs/
                else
                        if ((verbose)); then
                                echo "\nStarted genome assembly\n"
                        fi
			
			ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa mkdir ${outputFolder}gwa_output
			if ((kmersize)); then
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${unicyclerPath} -1 ${pathToInputFiles}/gwa_1.fq.gz -2 ${pathToInputFiles}/gwa_2.fq.gz --kmer_count $kmersize -t ${threads} -o ${outputFolder}gwa_output
			else
				ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw | xargs -I gwa ${unicyclerPath} -1 ${pathToInputFiles}/gwa_1.fq.gz -2 ${pathToInputFiles}/gwa_2.fq.gz -t ${threads} -o ${outputFolder}gwa_output
			fi
			rm -r ${pathToInputFiles}trimmed_reads
			for v in `ls ${pathToInputFiles} | grep _1.fq.gz | xargs -I gw basename -s _1.fq.gz gw`
			do
        			cp ${outputFolder}${v}_output/assembly.fasta ${outputFolder}assembled_outputs/${v}_assembled.fasta
        			rm -r ${outputFolder}${v}_output
			done
                fi
        else
                echo "Wrong input option for genome assembler. Type -h option for help"
                exit 1
        fi
}


main() {
        # Function that defines the order in which functions will be called

        get_input "$@"
        check_files
        if ((qualityControl)); then
                quality_control
        fi
        genome_assembly
}

# Calling the main function
main "$@"
