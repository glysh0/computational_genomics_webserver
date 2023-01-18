#!/bin/env python3
import argparse
import subprocess
import os

def chewBBACA(input_genomes,output_dir, cpu):
    '''Creates a wgMLST and cgMLST schema as well as a Newick tree of the cgMLST
       Input:
            input_genomes = directory where complete or draft genomes are located
            output_dir = directory for output
            cpu = Number of cpus
    '''
    os.chdir(output_dir)
    # Get training file
    if not os.path.isfile(output_dir+"prodigal_training_files/Escherichia_coli.trn"):
        subprocess.run("git clone https://github.com/mickaelsilva/prodigal_training_files", shell = True)
    import pdb; pdb.set_trace()
    # wgMLST Schema Creation
    wg = "chewBBACA.py CreateSchema -i " + input_genomes + " --cpu " + str(cpu) + " -o "+ output_dir + "Schema/" + " --ptf " "Escherichia_coli.trn"
    #print(wg.split())
    subprocess.run(wg.split())

    # Allele Calling
    allele = "chewBBACA.py AlleleCall -i " + input_genomes + " -g Schema/ -o "+ output_dir +"results_allele --cpu "+str(cpu)+" --ptf prodigal_training_files/Escherichia_coli.trn"
    #print(allele.split())
    subprocess.run(allele.split())

    # Define cgMLST schema
    date = subprocess.check_output("ls -t results_allele/ | head -n1", shell = True)
    date = str(date,'utf-8').rstrip()
    #print(date)
    cg = "chewBBACA.py ExtractCgMLST -i results_allele/"+ date +"/results_alleles.tsv -o " + output_dir +"cgMLST/ -r " + output_dir + "results_allele/" + date + "/RepeatedLoci.txt -p 0.95"
    #print(cg)
    subprocess.run(cg.split())

    # Create Newick Tree
    grapetree = "grapetree --profile " + output_dir + "cgMLST/cgMLST.tsv > cgMLST.tree"
    #print(grapetree)
    subprocess.run(grapetree, shell = True)

    return None

def collect_assembled_genomes(pathToInputFiles):
    #pathToInputFiles = "/home/projects/group-a/Team1-GenomeAssembly/assembled_output/"
    input_files = os.listdir(pathToInputFiles)
    base_names = []
    cleaned_files = []
    for file in input_files:
        if file.endswith(".fasta"):
            base_names.append(file)
            cleaned_files.append(pathToInputFiles+file)
    return(base_names, cleaned_files)

def MUMmer(prefix, reference_file, query_file):
    nucmer_command = ["nucmer", "-p", prefix, reference_file, query_file]
    delta_file = prefix+".delta"
    dnadiff_command = ["dnadiff", "-d", delta_file]
    subprocess.call(nucmer_command)
    subprocess.call(dnadiff_command)
    return(delta_file)

def kSNP(pathToInputDirectory,outDir,k,tree_type):

    '''
        - pathToInputDirectory: contains the path to the input directory where the fasta files are present
        - outDir: contains the path to where the output files are to be stored
        - k: User can select the k-mer size; when the k-mer size is not specified, the function calculates the k-mer size.
        - tree_type: should contain P, ML or NJ (only one). P for Parsimony Trees, ML for Maximum Likelihood Trees, and NJ for Neighbor-Joining Trees
    '''

    #creating a temp directory
    make = "mkdir snp_temp"
    snp_temp_path = os.path.abspath('snp_temp')
    subprocess.call(make.split())

    #getting the input directory name without its absolute path
    input_Directory = pathToInputDirectory.split("/")[-1]
    dir_path = pathToInputDirectory.strip(input_Directory)

    #storing the current working directory
    cwd = os.getcwd()

    #changing to the input directory to run MakeKSNP3infile
    os.chdir(dir_path)

    ################ creating an input list of fasta files with the respective genome ids ################
    input_list = "MakeKSNP3infile "+input_Directory+" "+snp_temp_path+"/kSNP_input.txt A"
    subprocess.call(input_list.split())

    #changing to the current working directory
    os.chdir(cwd)

    if k == 1: #default value, k-mer size for the dataset will be determined
        ############# running Kchooser to determine the k-mer size if not specified #########

        ######### creating a single multi-fasta file from the input files for kChooser #########
        makeFasta_command = "MakeFasta "+snp_temp_path+"/kSNP_input.txt "+snp_temp_path+"/kchooser_input.fasta"
        subprocess.call(makeFasta_command.split())

        ''' Note: The Kchooser program that comes along with the kSNP tool sometimes throws errors on certain systems.
        If the error is similar to: "system jellyfish dump output_0 -c > jellyout.txt failed: Inappropriate ioctl for device at Kchooser.pl line 242."
        Please download the corrected Kchooser script from: https://sourceforge.net/p/ksnp/discussion/general/thread/8f02db67/f844'''

        kchooser_command = "Kchooser "+snp_temp_path+"/kchooser_input.fasta"
        subprocess.call(kchooser_command.split())

        fhand = open('Kchooser.report','r')

        for line in fhand:
            line.rstrip()
            if line.startswith('The optimum value of K is'):
                line = line.split()
                k = line[6].strip(".")

        fhand.close()

        ########## running kSNP and visualizing the phylogenetic tree ################

    if tree_type == 'P':

        kSNP_command = "kSNP3 -in "+snp_temp_path+"/kSNP_input.txt -outdir "+outDir+" -k "+str(k)

        figtree_command = "figtree -graphic PDF "+outDir+"/tree.parsimony.tre "+outDir+"/kSNP_run_parsimony.pdf"

    if tree_type == 'ML':

        kSNP_command = "kSNP3 -in "+snp_temp_path+"/kSNP_input.txt -outdir "+outDir+" -k "+str(k)+" -ML"

        figtree_command = "figtree -graphic PDF "+outDir+"/tree.ML.tre "+outDir+"/kSNP_run_ML.pdf"

    if tree_type == 'NJ':

        kSNP_command = "kSNP3 -in "+snp_temp_path+"/kSNP_input.txt -outdir "+outDir+" -k "+str(k)+" -NJ"

        figtree_command = "figtree -graphic PDF "+outDir+"/tree.NJ.tre "+outDir+"/kSNP_run_NJ.pdf"

    subprocess.call(kSNP_command.split())
    subprocess.call(figtree_command.split())

    #deleting the temporary directory
    delete = "rm -r snp_temp"
    subprocess.call(delete.split())

    return None

def main():

    #get options using os or argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="directory for input genomes", required = True)
    parser.add_argument("-o", "--output", help="name of output directory", required = True)
    parser.add_argument("-c", "--cpu", help="Number of cpus to use", default = 6)
    parser.add_argument("-t", "--tool", help="Tool of choice: MUMmer (m), chewBBACA (c), kSNP (k), or all (a)", choices = ['m', 'c', 'k', 'a'], required = True)
    parser.add_argument("-k","--kmersize",help="Specify K-mer Size for kSNP",default=1)
    parser.add_argument("-tr","--treetype",help="Phylogenetic Tree of Choice for kSNP output", choices = ['P','ML','NJ'], default='P')
    args = parser.parse_args()

    output = args.output
    if output[-1] != "/":
        output += "/"

    #call MUMmer
    names, mummer_inputs = collect_assembled_genomes(args.input)
    reference = mummer_inputs[0]

    for i in range(0, len(mummer_inputs)):
        prefix = args.input + names[i]
        MUMmer(prefix, reference, mummer_inputs[i])

    #call chewBBACA
    if args.tool == 'c' or args.tool == 'a':
        chewBBACA(args.input, output, args.cpu)

    #call kSNP

    k = args.kmersize
    tr = args.treetype

    kSNP(args.input,args.output,k,tr)

    #analysis or visualization -- if we're including it here

if __name__ == "__main__":
    main()
