import os
import csv
import shutil
import logging
import argparse
import datetime
import subprocess
import pandas as pd
from datetime import datetime
from Bio import SeqIO

'''
The Specie Detection check if the raw data species (fastq) we got from WGS
It run kraken2 to find the species and create a report.
It's help us to compare them to the file name given after the WGS (ECXXX for Ecoli
NMXXX for n. meningitidis etc)
To discriminate EColi ans Shigella, we need to use another workflow: de novo assembly and 
use blast to check the presence of lacY, a gene present in all EColi but not in Shigella
'''

def if_no_create_dir(d):
    try:
        os.makedirs(d, exist_ok=True)
    except Exception as e:
        exit(f'Cannot create the folder {d}\n{e}')



def gather_reads(arr):
    paired_files = {}
    extensions = ['.fastq', '.fq', '.fasq', '.fq.gz', '.fastq.gz']

    for file_name in arr:
        if any(file_name.endswith(ext) for ext in extensions):
            if '_R1' in file_name:
                sample_name = file_name.split('_R1')[0]
                paired_files.setdefault(sample_name, {}).update({'R1': file_name})
            elif '_R2' in file_name:
                sample_name = file_name.split('_R2')[0]
                paired_files.setdefault(sample_name, {}).update({'R2': file_name})

    return paired_files

def get_specie(file):
    cmd = "awk -F'\t' '$4 == \"S\" {print $6; exit}' " + file

    output = subprocess.check_output(cmd, shell=True, text=True)
    values = output.strip().split()

    if len(values) == 2:
        genus, specie = values
        return genus + ' ' + specie
    elif len(values) == 1:
        return values[0]
    else:
        return ''

def create_specie_detection_report(inp, out):
    specie_detected = get_specie(inp) 

    if not os.path.exists(out):
        try:
            os.makedirs(out)
        except:
            pass

    now = datetime.now()
    formatted_date = now.strftime("%Y%m%d")

    report_csv = f'{my_output}/specie_detection_{formatted_date}.csv' 

    if not os.path.exists(report_csv):
        with open(report_csv, 'w') as f:
            f.write(f'sample\tSpecie Detected\n')
            f.write(f'{sample}\t{specie_detected}\n')
    else:
        with open(report_csv,"a") as f:
            f.write(f'{sample}\t{specie_detected}\n')

def fetch_ecoli_and_shigella(inp, x):
    matching = []
    others   = []
    for file_name in os.listdir(inp):
        if any(substring in file_name for substring in x):
            matching.append(file_name)
        else:
            others.append(file_name)

    return matching, others

def trimming(r1, r2, inp, out, cpus=8):
    cmd = f"""
    trimmomatic PE -phred33 -threads {cpus} \
    {inp}/{r1} {inp}/{r2} \
    {out}/{r1} {out}/trim_unpaired_{r1} \
    {out}/{r2} {out}/trim_unpaired_{r2} \
    HEADCROP:20 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 CROP:265 MINLEN:50;
    rm -f {out}/trim_unpaired_*;
    """
    os.system(cmd)

def metagenomics(r1, r2, output, sample, db, cpus=10):
    try:
        cmd = f"""
        kraken2 --use-names --db {db} \
        --threads {cpus} \
        --report {output}/{sample}.report \
        --paired {r1} {r2} \
        --output - 2>&1
        """
        os.system(cmd)
    except:
        with open(f'{sample}.error', 'w') as f:
            f.write(f'Cannot get the specie\n')

def assembly(r1, r2, output, sample, cpus=10):
    cmd =  f'spades.py -1 {r1} -2 {r2} -t {cpus} -o {output} 2>&1;'
    cmd += f'mv {output}/scaffolds.fasta {output}/${sample}.fasta'

    os.system(cmd)

def modify_fasta_header(fasta_file, tmp_output):
    name = os.path.splitext(os.path.basename(fasta_file))[0]
    cmd = f"awk -v name='{name}' '{{if ($0 ~ /^>/) print \">\" name; else print $1;}}' {fasta_file} >> {tmp_output}/{name}.fasta"
    
    os.system(cmd)

def combine_fastas(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    tmp_dir = 'tmp'
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    for file in os.listdir(input_dir):
        modify_fasta_header(f'{input_dir}/{file}', tmp_dir)
        os.system(f'cat {tmp_dir}/{file} >> {seq_combined_file}')
    
        
def create_db(file, output_dir):
    os.system(f'makeblastdb -in {file} -dbtype nucl -out {output_dir}/db')

def blast(query, db, out):
    os.system(f'blastn -query {query} -db {db} -out {out} -outfmt 6')

def filter_blast_results(
        input_dir, output_csv, min_length=1200, 
        min_identity=90.0, max_evalue=1e-5):
    filtered_results = []
    
    header = ['subject', SAMPLE_NAME_COL, 'perc_identity', 'alignment_len', 'mismatches', 'gaps', 
              'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', GENE_NAME_COL]
    
    # Iterate over all files in the input directory
    for file_name in os.listdir(input_dir):
        gene = file_name.split('.')[0]
        # Only process .txt files (assuming BLAST results are stored as txt)
        if file_name.endswith('.txt'):
            file_path = os.path.join(input_dir, file_name)
            
            with open(file_path, 'r') as file:
                for line in file:
                    # Skip empty lines or headers
                    if line.strip() == "" or line.startswith("#"):
                        continue
                    
                    # Split the line by tabs
                    columns = line.strip().split('\t')
                    
                    # Ensure there are enough columns to parse
                    if len(columns) < 12:
                        continue
                    
                    # Parse fields
                    subject, query = columns[0], columns[1]
                    perc_identity = float(columns[2])
                    alignment_len = int(columns[3])
                    mismatches = columns[4]
                    gaps = columns[5]
                    q_start = columns[6]
                    q_end = columns[7]
                    s_start = columns[8]
                    s_end = columns[9]
                    e_value = float(columns[10])
                    bit_score = float(columns[11])
                    
                    # Apply filtering criteria
                    if alignment_len >= min_length and perc_identity >= min_identity and e_value <= max_evalue:
                        # Append the filtered row with the source file name
                        filtered_results.append([
                            subject, query, perc_identity, 
                            alignment_len, mismatches, gaps, 
                            q_start, q_end, s_start, s_end, e_value, 
                            bit_score, gene
                        ])
    
    with open(output_csv, 'w', newline='') as out_file:
        csv_writer = csv.writer(out_file)
        csv_writer.writerow(header) 
        csv_writer.writerows(filtered_results)  
######################################################
#####################################################
# Args

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Here your genomes dir")
parser.add_argument("--output", default='results', required=False, help="Here your destination directory")
args = parser.parse_args()


######################################################
#####################################################
# Variables

DB                     = '/home/user1/Desktop/analysis/genome_analysis/dev/bin/minikraken2_v2_8GB_201904_UPDATE'
GENES_DIR              = 'genes'
TRIMMED_DIR            = 'trimmed'
ASSEMBLIES_DIR         = 'assemblies'
SEQ_COMBINED_DIR       = 'sequences_combined'
DB_DIR                 = 'db'
BLAST_RESULTS_DIR      = 'blast_results'
SPECIE_DETECTION_DIR   = 'specie_detection'

LOGS_FILE              = 'logs.txt'
EC_SG_REPORT_FILE      = 'ec_sg_report.csv'
EC_SG_SUM_REPORT_FILE  = 'ec_sg_sum_report.csv'
DETECTION_REPORT_FILE  = 'detection_report.csv'

GENE_NAME_COL          = 'gene'
SAMPLE_NAME_COL        = 'sample'

now = datetime.now()
formatted_date = now.strftime("%Y%m%d_%H%M%S")


genomes             = args.input
wd                  = os.path.dirname(os.path.realpath(__file__))
genes_dir           = os.path.join(wd, GENES_DIR)
output_dir          = os.path.join(wd, args.output)
assemblies_dir      = os.path.join(wd, ASSEMBLIES_DIR)

specie_detection_dir  = os.path.join(args.output, f'{SPECIE_DETECTION_DIR}_{formatted_date}')
db_dir                = os.path.join(specie_detection_dir, DB_DIR)
blast_dir             = os.path.join(specie_detection_dir, BLAST_RESULTS_DIR)
seq_combined_dir      = os.path.join(specie_detection_dir, SEQ_COMBINED_DIR)

db_name             = os.path.join(db_dir, 'db')

logs_file            = os.path.join(specie_detection_dir, LOGS_FILE)
seq_combined_file   = os.path.join(seq_combined_dir, 'seq_combined_file')
blast_report_file   = os.path.join(blast_dir, EC_SG_REPORT_FILE)
blast_sum_report_file   = os.path.join(blast_dir, EC_SG_SUM_REPORT_FILE)
detection_report_file  = os.path.join(specie_detection_dir, DETECTION_REPORT_FILE)
ec_and_sh           = ['EC', 'F-EC', 'SG']


logging.basicConfig(filename=logs_file, level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

#############################################################
############################################################
# Start

if not os.path.exists(specie_detection_dir):
    os.makedirs(specie_detection_dir, exist_ok=True)

arr_ecoli_shigella_reads, arr_others_reads = fetch_ecoli_and_shigella(genomes, ec_and_sh)

ec_sh_species   = gather_reads(arr_ecoli_shigella_reads)
others_species  = gather_reads(arr_others_reads)


if ec_sh_species:
    combine_fastas(assemblies_dir, seq_combined_dir)
    
    if_no_create_dir(db_dir)
    if_no_create_dir(blast_dir)
    
    create_db(seq_combined_file, db_dir)

    for d in os.listdir(genes_dir):
        for f in os.listdir(f'{genes_dir}/{d}'):
            blast(
                query=os.path.join(genes_dir, d, f),
                db=f'{db_dir}/db',
                out=os.path.join(blast_dir, f'{f.split(".")[0]}.txt')
            )

    if os.path.exists(blast_dir):
        filter_blast_results(blast_dir, blast_report_file)
    
    # Check the specie - gene match
    if os.path.exists(blast_report_file):
        df = pd.read_csv(blast_report_file)
    
    values = []
    for d in os.listdir(genes_dir):
        for f in os.listdir(f'{genes_dir}/{d}'):
            gene_name = f.split('.')[0]
            if df[GENE_NAME_COL].isin([gene_name]).any():
                sample = df.loc[df[GENE_NAME_COL] == gene_name, SAMPLE_NAME_COL].values[0]
                values.append([sample , d])     # the specie name is the folder name (d)

    if values:
        header = ['sample', 'specie detected']
        with open(blast_sum_report_file, 'w', newline='') as out_file:
            csv_writer = csv.writer(out_file)
            csv_writer.writerow(header) 
            csv_writer.writerows(values) 
    else:
        logging.error(f"[EC_SG Detection] Error occurred while fetching report values\n")


#     for specie in arr_ecoli_shigella_r1:
#         sample = specie.split('_')[0]
#         ext = specie.split('.')[1]

#         inp = os.path.join(output_dir, sample, TRIMMED_DIR)
#         out = os.path.join(output_dir, sample, _ES)

#         assembly(
#             r1=f'{inp}/{sample}_R1.{ext}',
#             r2=f'{inp}/{sample}_R2.{ext}',
#             sample=sample,
#             output=out,
#         )   

# if others_species:
#     for sample, files in others_species.items():
#         r1 = files['R1']
#         r2 = files['R2']

#         # inp = os.path.join(output_dir, sample, TRIMMED_DIR)
#         inp = genomes
        
#         metagenomics(
#             r1=f'{inp}/{r1}',
#             r2=f'{inp}/{r2}',
#             sample=sample,
#             output=specie_detection_dir,
#             db=DB
#         )




