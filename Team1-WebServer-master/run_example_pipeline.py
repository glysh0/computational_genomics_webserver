#!/projects/VirtualHost/predicta/html/miniconda3/bin/python3.7

import os
import glob
import subprocess
import tempfile
import sys
import argparse
import yagmail
import shutil

class Pipeline:
    """
    Runs pipeline either from beginning or intermediary stage
    """
    def __init__(self, tmp_folder, input_path, epidata_path,
                 assembly_parameters,
                 gene_prediction_parameters,
                 functional_annotation_parameters,
                 comparative_genomics_parameters):
        """
        Initialize class and run demanded steps of pipeline
        tmp_folder: str, path to temporary folder
        assemble: boolean, Perform assembly of contigs
        predict_genes: boolean, perform gene prediction
        functional_annotation: boolean, perform functional annotation
        comparative_genomics: boolean, perform comparative genomics
        """
        self.tmp_folder = tmp_folder
        self.input_path = input_path
        self.epidata_path = epidata_path
        self.assembly_parameters = assembly_parameters
        self.gene_prediction_parameters = gene_prediction_parameters
        self.functional_annotation_parameters = functional_annotation_parameters
        self.comparative_genomics_parameters = comparative_genomics_parameters
        # initialize results
        self.assembly = None
        self.gene_prediction = None
        self.functional_annotation = None
        self.comparative_genomics = None
        # run assembling
        if self.assembly_parameters['assemble']:
           self.assembly = self.run_genome_assembly()

        # run gene prediction
        if self.gene_prediction_parameters['predict_genes']:
            if self.assembly is None:
                input_path = self.input_path
            else:
                input_path = self.tmp_folder + 'assembled_outputs/'
            self.gene_prediction = self.run_gene_prediction(input_path)
        # run functional annotation
        if self.functional_annotation_parameters['functional_annotation']:
            if self.gene_prediction is None:
                #TODO: remove amino acids
                input_path = self.input_path + 'amino_acids/'
            else:
                input_path = self.tmp_folder + '/amino_acids/'
            self.functional_annotation = self.run_functional_annotation(input_path)

        # run comparative_genomics:
        if self.comparative_genomics_parameters['comparative_genomics']:
            if self.functional_annotation is None:
                input_path = self.input_path
            else:
                input_path = self.tmp_folder + 'assembled_outputs/'
            self.comparative_genomics = self.run_comparative_genomics(input_path)


    def run_genome_assembly(self):
	# output2 = subprocess.getoutput('/home/dkesar3/Team1-WebServer/scripts/quast.sh -t 8 -q /home/projects/group-a/bin/quast/quast.py -p test_data -u /home/dkesar3/Team1-WebServer/unicycler_contigs -o /home/dkesar3/Team1-WebServer/assembled_contigs/ -v')
        spades = self.assembly_parameters['spades']
        trimming = self.assembly_parameters['trimming']
        kmer_size = self.assembly_parameters['kmer_size']        
        #self.output_path_assembly =f'{self.tmp_folder}'
        log_file = open(f'{self.tmp_folder}/genomeAssemblyLog.txt','w+')
        options = []
        #TODO: Add trimming option
        if not spades:
            options.append('-g')
            options.append('a')
        else:
            options.append('-g')
            options.append('s')
        if trimming:
            options.append('-q')
            options.append('1')
        if not kmer_size is None:
            options.append('-k')
            options.append(str(kmer_size))
        # run unicycler
        output = subprocess.check_output([f"{cwd}scripts/run_unicycler.sh", "-t", "8", "-p", self.input_path, "-o", self.tmp_folder, "-v", '1'] + options)
        quast_output = subprocess.check_output([f"{cwd}scripts/run_quast.sh", "-t", "8", "-p", self.input_path, "-o", self.tmp_folder, "-q", "quast.py", "-v"])

        log_file.write(str(output))
        log_file.close()
        return output

    def run_gene_prediction(self, input_path):
        options = []
        if self.gene_prediction_parameters['gms2']:
            options.append('-g')
        elif self.gene_prediction_parameters['trnascan']:
            options.append('-t')
        elif self.gene_prediction_parameters['rnammer']:
            options.append('-r')
        
        log_file = open(f'{self.tmp_folder}/genePredictionLog.txt','w+')
        # run gene prediction on all fasta files
        assembled_files = [f for f in glob.glob(f'{input_path}/*.fasta') if os.path.isfile(f)]
        for l in assembled_files:
            input_dir = [f"{cwd}scripts/run_dfast.sh", "-i",\
                         l, "-o", self.tmp_folder +'/', "-v" ] + options
            output = subprocess.check_output(input_dir)
            log_file.write(str(output))
        log_file.close()
        return output

    def run_functional_annotation(self, input_path):
        log_file = open(f'{self.tmp_folder}/funtionalAnnotationLog.txt', 'w+')
        options = []
        options.append('-u')
        options.append('usearch11.0.667_i86linux32')
        if self.functional_annotation_parameters['eggnog']:
            options.append("-e")
            options.append("/projects/VirtualHost/predicta/html/Team1-WebServer/tools/eggnog-mapper/emapper.py")
        if self.functional_annotation_parameters['signalp']:
            options.append("-p")
            options.append("signalp")
        if self.functional_annotation_parameters['tmhmm']:
             options.append("-t")
             options.append("tmhmm")
        if self.functional_annotation_parameters["interpro"]:
             options.append("-s")
             options.append("interproscan.sh")
        if self.functional_annotation_parameters['deeparg']:
             options.append("-d")
             options.append("/projects/VirtualHost/predicta/html/Team1-WebServer/tools/deeparg-ss/deepARG.py")

        cmd  = [f"{cwd}scripts/functional_annotation.py", "-f", input_path, "-o", self.tmp_folder] + options
        output = subprocess.check_output(cmd)
        log_file.write(str(output))
        log_file.close()
        return output

    def run_comparative_genomics(self, input_path):
        log_file = open(f'{self.tmp_folder}/comparativeGenomicsLog.txt', 'w+')
        cmd = [f"{cwd}scripts/cg_pipeline.py", '-i', input_path, '-o', self.tmp_folder, '-c', '8', '-t', self.comparative_genomics_parameters['tools']]
        output = subprocess.check_output(cmd)
        log_file.write(str(output))
        log_file.close()
        return output


class Results:
    """
    Class that extracts results and prepares email send to user if provided
    """
    def __init__(self, tmp_folder, user_email):
        self.tmp_folder = tmp_folder
        self.job_id = tmp_folder.split('/')[-1]
        self.user_email = user_email
        self.files_to_send = []
        #import pdb; pdb.set_trace()
        if os.path.isfile(f'{self.tmp_folder}/assembly_statistics.txt'):
            self.files_to_send.append(f'{self.tmp_folder}/genomeAssemblyLog.txt')
        if len(glob.glob(f'{self.tmp_folder}/*.gff')) >= 1:
             predictions = glob.glob(f'{self.tmp_folder}/*.gff')
             self.extract_gene_prediction_statistics(predictions)
             self.files_to_send.append(f"{self.tmp_folder}/gene_prediction_statistics.txt")
        #self.zipped_files = zipfile.ZipFile(f"{self.tmp_folder}/log_files.zip", 'a')
        #for f in files_to_zip:
        #    self.zipped_files.write(f, compress_type = zipfile.ZIP_DEFLATED)
        self.send_email()
        #TODO:extensions of fa and cg files?
        #self.functional_annoation = os.listdir(f'{self.tmp_folder}/*.?')
        #self.comparative_genomics = os.listdir(f'{self.tmp_folder}/*.?')
    #TODO what do we want to sent? which files --> will need to compress them


    def extract_gene_prediction_statistics(self, list_of_gff):
        with open(self.tmp_folder + '/gene_prediction_statistics.txt' , 'w+') as out:
            out.write("Sample #CDS #tRNA #rRNA #CRISPR\n")
            for p in list_of_gff:
                statistics = subprocess.check_output([f'{cwd}extract_statistics_gene_prediction.sh', p])
                out.write(str(statistics.strip())[2:-1])
        out.close()


    def send_email(self):
        body = f"Results from ECHO for job {self.job_id}"
        # use oauth2 for authorization requires oauth2.json file to be in same folder as cwd
        yag = yagmail.SMTP("echowebserver@gmail.com", oauth2_file=f'{cwd}oauth2.json')
        yag.send(
        to=self.user_email,
        subject="ECHO results",
        contents=body, 
        attachments=self.files_to_send,
        )


def start_to_end(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assemble', help='Assemble contigs, requires pair-end-reads in fastq format, default=0',
                        action='store_true', default=False, dest='assemble', required=False)
    parser.add_argument('--spades', help='Run SPADES instead of Unicycler, default=0', default=False, action='store_true',
                        required=False, dest='spades')
    parser.add_argument('--trim', help='Perform trimming, default=0', default=False, action='store_true', required=False, dest='trim')
    parser.add_argument('--kmer_size', help='Kmer-size used by unicycler', default=None, required=False, dest='kmer_size')
    parser.add_argument('-p', '--predict_genes', help='Predict genes in assembled contigs, requires either files in FASTA format or assembly must be run',
                        action='store_true', default=False, dest='predict_genes', required=False)
    parser.add_argument('--gms2', help='Use GeneMarkS2 for CDS prediction instead of Prodigal', default=False,
                        action='store_true', required=False, dest='gms2')
    parser.add_argument('--trnascan', help='Use tRNAscan-SE instead for tRNA prediciton instead of Aragorn',
                        default=False, action='store_true', dest='trnascan')
    parser.add_argument('--rnammer', help='Run RNAmmer for rRNA prediciton instead of barrnap',
                        default=False, action='store_true', dest='rnammer')
    parser.add_argument('-f', '--functional_annoation', help='Perform functional annotation, ',
                        action='store_true', default=False, dest='functional_annotation', required=False)
    #parser.add_argument('--usearch', help='Run Usearch', default=False, action='store_true', dest='usearch', required=False)
    parser.add_argument('--eggnog', help='Run eggnog', default=False, action='store_true', dest='eggnog', required=False)
    parser.add_argument('--tmhmm', help='Run tmhmm', default=False, action='store_true', dest='tmhmm', required=False)
    parser.add_argument('--signalp', help='Run signalp', default=False, action='store_true', dest='signalp', required=False)
    parser.add_argument('--deeparg', help='Run DeepARG', default=False, action='store_true', dest='deeparg', required=False)
    parser.add_argument('--interpro', help='Run interproscan', default=False, action='store_true', dest='interpro', required=False)

    parser.add_argument('-c', '--comparative_genomics', help='Perform comparative genomics, ',
                        action='store_true', default=False, dest='comparative_genomics', required=False)
    parser.add_argument('--cg_tools', help='tool of choice: MUMmer (m), chewBBACA (c), kSNP3 (k), all (a), default=a',
                        default='a', choices = ['m', 'c', 'k', 'a'], required=False, dest='cg_tools')
    parser.add_argument('-e', '--email', help="Email address of user to which results will be send",
                        default=None, required=False, dest='email')
    parser.add_argument('-i', '--input', help='Path to input files', required=True, dest='input')
    parser.add_argument('--epidata', help='Path to epidata', required=False, default=None, dest='epidata')
    args = parser.parse_args()
    args = vars(args)
    input_path = args['input']

    assembly_parameters = {'assemble': args['assemble'],
                           'spades': args['spades'],
                           'trimming': args['trim'],
                           'kmer_size': args['kmer_size']}

    gene_prediction_parameters = {'predict_genes': args['predict_genes'],
                                  'gms2': args['gms2'],
                                  'trnascan': args['trnascan'],
                                  'rnammer': args['rnammer']}
    
    functional_annotation_parameters = {'functional_annotation': args['functional_annotation'],
                                        'eggnog': args['eggnog'],
                                        'tmhmm': args['tmhmm'],
                                        'signalp': args['signalp'],
                                        'deeparg': args['deeparg'],
                                        'interpro': args['interpro']}

    comparative_genomics_parameters = {'comparative_genomics': args['comparative_genomics'],
                                      'tools': args['cg_tools']}

    user_email = args['email']
    epidata = args['epidata']
    
    global cwd
    cwd = os.getcwd() + '/'
    # create tmp dir in cwd
    current_tmp_dir = tempfile.mkdtemp(prefix=cwd + 'analysis/')
    #current_tmp_dir = '/projects/VirtualHost/predicta/html/Team1-WebServer/analysis/xgn_dsdc/'
    #subprocess.call(["python", "./preinstall.py"])
    pipeline = Pipeline(current_tmp_dir, input_path, epidata,
                        assembly_parameters,
                        gene_prediction_parameters,
                        functional_annotation_parameters,
                        comparative_genomics_parameters)
    if user_email is not None:
        Results(current_tmp_dir, user_email)
    # done with everything --> clean up
    #shutil.rmtree(current_tmp_dir)

if __name__ == "__main__":
    start_to_end(sys.argv[1:])
