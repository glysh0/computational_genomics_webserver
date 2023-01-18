#! /usr/bin/env python
from flask import Flask, render_template, request, jsonify, flash
import os
import subprocess
import sys

app = Flask(__name__,
            static_url_path='', 
            static_folder='templates/'
            )
app.config['SECRET_KEY'] = 'WEBSERVERSECRETKEY'

@app.route('/General_Pipeline', methods=['GET','POST'])
def General_Pipeline():
    if request.method == "POST":
        #print (request.form)
        #print (request.files)
        email = request.form.get('emailaddress')
        file = request.files.get('file1')
        metadata = request.files.get('metadata')
        filename = file.filename 
        #metadataname = metadata.filename
        upload_dir = os.path.join(os.getcwd(), 'upload/')
        if not os.path.exists(upload_dir):
            os.makedirs(upload_dir)
        filePath = os.path.join(upload_dir, filename)
        #metadataPath = os.path.join(upload_dir, metadataname)
        file.save(filePath)
        #metadata.save(metadataPath)
        #run example_pipline
        return_info = subprocess.Popen([sys.executable, "./run_example_pipeline.py", "-a","--trim","-p","-f","-c","--cg_tools", "a", "-e", email, "-i", filePath], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            next_line = return_info.stdout.readline()
            return_line = next_line.decode("utf-8", "ignore")
            if return_line == '' and return_info.poll() != None:
                break
            if return_line:
                print(return_line)
        
        returncode = return_info.wait()
        if returncode:
            print (subprocess.CalledProcessError(returncode, return_info))

    return render_template('analyze.html')

@app.route('/home')
def home():
    return render_template('home.html')
@app.route('/team')
def team():
    return render_template('team.html')

@app.route('/class_results')
def class_results():
    return render_template('class_results.html')


if __name__=="__main__":
    #app.run(host="biogenome2020.biosci.gatech.edu",debug=True)
    app.run(host="predict2020t1.biosci.gatech.edu",debug=True)

