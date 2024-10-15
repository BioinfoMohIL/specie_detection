import os
import re
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

def get_info(step, msg):
    t = f'\n{step}          {msg}'
    print(t)
    logging.info(t)

def error_handler(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in {func.__name__}: {e}", exc_info=True)
            exit()
    return wrapper

def create_dir(d):
    try:
        os.makedirs(d, exist_ok=True)
        logging.info(f"[create_dir] Successfully created directory: {d}\n")
    except Exception as e:
        msg = f'Cannot create the folder {d}\n{e}'
        logging.error(f"[create_dir] Error occurred while creating directory {d}: {e}\n")
        exit(msg)

def remove_dir(dir_path):
    try:
        shutil.rmtree(dir_path)
        get_info('[Clean up]', f"Successfully removed the directory: {dir_path}")
    except OSError as e:
        get_info('[Clean up]', f"Error: {e.strerror}")

def unzipReads(file):
    cmd = f'gunzip -f {file}'
    os.system(cmd)

def combineFiles(dir):
    cmd ='cd ' + dir + ';'
    cmd+= 'ls ' + dir + ' | egrep "_L00" | '
    cmd+= 'awk \'BEGIN{FS=OFS="_"}{ print $1 } \' | sort | uniq | '
    cmd+= 'while    read    line;do cat ' + dir + '/${line}_*_*_R1*.fastq.gz > ${line}_R1.fastq.gz && cat ' + dir + '/${line}_*_*_R2*.fastq.gz > ${line}_R2.fastq.gz; done;'
    cmd+= 'cd ..'
    os.system(cmd)

def renameFile(dir):
    for file in os.listdir(dir):
        sample_name = file.split('_')[0]

        if '_R1' in file:
            sample_name += '_R1'
        elif '_R2' in file:
            sample_name += '_R2'

        if '.gz' in file:
            extension = 'fastq.gz'
        else:
            extension = 'fastq'

        name = f"{sample_name}.{extension}"
        
        if(file == name):
           continue

        cmd = f'mv {os.path.join(dir, file)} {os.path.join(dir, name)}'
        os.system(cmd) 

@error_handler
def setup(input_dir):
    filenames = os.listdir(input_dir)

    pattern = r"^(.*)_R[12](?:_\d+)?\.(fastq|fq)(\.gz)?$"

    for filename in filenames:
        match = re.match(pattern, filename)
        if match:
            reads_name = match.group(1) or match.group(2)
        else:
            exit(f"Filename: {filename}, Read unknow")

    isMoreThanTwoReads = [s for s in os.listdir(input_dir) if "_L002" in s or "_L003" in s or "_L004" in s]

    if(len(isMoreThanTwoReads) > 0):
        try:
            combineFiles(input_dir);
        except Exception as e:
            exit(f"Cannot rename reads. Please check your reads\n{e}")
   
    renameFile(input_dir)
    
    for file in os.listdir(input_dir):
        if '.gz' in file:
            try:
                unzipReads(f'{input_dir}/*.gz')
                break
            except Exception as e:
                exit(f"Cannot unzip {file}. Please check your reads\n{e}")  

@error_handler
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

@error_handler
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
    
@error_handler
def create_detection_report(input_file, output_file, sample):
    specie_detected = get_specie(input_file) 

    if not os.path.exists(input_file):
        with open(output_file, 'w') as f:
            f.write(f'{SAMPLE_NAME_COL},{SPECIE_DETECTED_COL}\n')
            f.write(f'{sample},{specie_detected}\n')
    else:
        with open(output_file,"a") as f:
            f.write(f'{sample},{specie_detected}\n')


def fetch_ecoli_and_shigella(inp, x):
    matching = []
    others   = []
    for file_name in os.listdir(inp):
        if any(substring in file_name for substring in x):
            matching.append(file_name)
        else:
            others.append(file_name)

    return matching, others

@error_handler
def trimming(r1, r2, inp, out, cpus=8):
    cmd = f"""
    trimmomatic PE -phred33 -threads {cpus} \
    {inp}/{r1} {inp}/{r2} \
    {out}/{r1} {out}/trim_unpaired_{r1} \
    {out}/{r2} {out}/trim_unpaired_{r2} \
    HEADCROP:20 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 CROP:265 MINLEN:50 > /dev/null 2>&1;
    rm -f {out}/trim_unpaired_*;
    """
    os.system(cmd)

@error_handler
def metagenomics(r1, r2, report, db, cpus=10):
    cmd = f"""
        kraken2 --use-names --db {db} \
        --threads {cpus} \
        --report {report} \
        --paired {r1} {r2} \
        --output - > /dev/null 2>&1
        """
    os.system(cmd)
   
@error_handler
def assembly(r1, r2, output, sample, cpus=30):
    cmd =  f'spades.py -1 {r1} -2 {r2} -t {cpus} -o {output} > /dev/null;'
    cmd += f'mv {output}/scaffolds.fasta {output}/{sample}.fasta'

    os.system(cmd)

@error_handler
def gather_assemblies(source, destination):
    shutil.move(source, destination)

def modify_fasta_header(fasta_file, tmp_output):
    name = os.path.splitext(os.path.basename(fasta_file))[0]
    cmd = f"awk -v name='{name}' '{{if ($0 ~ /^>/) print \">\" name; else print $1;}}' {fasta_file} >> {tmp_output}/{name}.fasta"
    
    os.system(cmd)

def combine_fastas(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    for file in os.listdir(input_dir):
        modify_fasta_header(f'{input_dir}/{file}', tmp_dir)
        os.system(f'cat {tmp_dir}/{file} >> {seq_combined_file}')

@error_handler         
def create_db(file, output_dir):
    os.system(f'makeblastdb -in {file} -dbtype nucl -out {output_dir}/db')

@error_handler
def blast(query, db, out):
    os.system(f'blastn -query {query} -db {db} -out {out} -outfmt 6 2>&1')

@error_handler
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

DB                     = 'path/to/minikraken2_v2_8GB_201904_UPDATE'
DB_DIR                 = 'db'
TMP_DIR                = 'tmp'
LOGS_DIR               = 'logs'
GENES_DIR              = 'genes'
EC_SG_DIR              = 'ec_sg'
SAMPLES_DIR            = 'samples'
TRIMMED_DIR            = 'trimmed'
ASSEMBLY_DIR           = 'assembly'
ASSEMBLIES_DIR         = 'assemblies'
SEQ_COMBINED_DIR       = 'sequences_combined'
BLAST_RESULTS_DIR      = 'blast_results'
SPECIE_DETECTION_DIR   = 'specie_detection'

LOGS_FILE              = 'logs.txt'
MATCHED_FILE           = 'matched.txt'
SEQ_COMBINED_FILE      = 'seq_combined_file'
EC_SG_REPORT_FILE      = 'ec_sg_report.csv'
EC_SG_SUM_REPORT_FILE  = 'ec_sg_sum_report.csv'
DETECTION_REPORT_FILE  = 'detection_report.csv'

GENE_NAME_COL          = 'gene'
SAMPLE_NAME_COL        = 'sample'
SPECIE_DETECTED_COL    = 'specie detected'
LACY                   = 'lacY' 
IPAH_1                 = 'ipaH_1'
IPAH_7                 = 'ipaH_7'

HEADER_SAMPLE          = 'sample'
HEADER_SPECIE          = 'specie detected'

now = datetime.now()
formatted_date = now.strftime("%Y%m%d_%H%M%S")

genomes             = args.input
wd                  = os.path.dirname(os.path.realpath(__file__))
genes_dir           = os.path.join(wd, GENES_DIR)
output_dir          = os.path.join(wd, args.output)
tmp_dir             = os.path.join(wd, TMP_DIR)

specie_detection_dir  = os.path.join(output_dir, f'{SPECIE_DETECTION_DIR}_{formatted_date}')
assemblies_dir        = os.path.join(specie_detection_dir, ASSEMBLIES_DIR)
logs_dir              = os.path.join(specie_detection_dir, LOGS_DIR)

ec_sg_dir             = os.path.join(specie_detection_dir, EC_SG_DIR)
samples_dir           = os.path.join(ec_sg_dir, SAMPLES_DIR)
db_dir                = os.path.join(ec_sg_dir, DB_DIR)
blast_dir             = os.path.join(ec_sg_dir, BLAST_RESULTS_DIR)
seq_combined_dir      = os.path.join(ec_sg_dir, SEQ_COMBINED_DIR)

db_name               = os.path.join(db_dir, 'db')

logs_file               = os.path.join(logs_dir, LOGS_FILE)
seq_combined_file       = os.path.join(seq_combined_dir, SEQ_COMBINED_FILE)
blast_report_file       = os.path.join(blast_dir, EC_SG_REPORT_FILE)
blast_sum_report_file   = os.path.join(blast_dir, EC_SG_SUM_REPORT_FILE)
detection_report_file   = os.path.join(specie_detection_dir, DETECTION_REPORT_FILE)
matched_file            = os.path.join(specie_detection_dir, MATCHED_FILE)

ec_and_sh           = ['EC', 'F-EC', 'SG']
ec_genes = [LACY]
sg_genes = [IPAH_1]
species_by_genes = {
    'lacY' : 'ecoli',
    'ipaH_1': 'shigella',
    'ipaH_7': 'shigella'
}

species_list = {
    'BP'    : 	'Bordetella pertussis',
    'EC'    : 	'Ecoli',
    'F-EC'  : 	'Ecoli',
    'SG'    : 	'Shigella',
    'HI'    : 	'Haemophilus influenzae',
    'LC'    : 	'Listeria monocytogenes',
    'LF'    : 	'Listeria monocytogenes',
    'LG'    : 	'Legionella pneumophila',
    'LW'    : 	'Legionella pneumophila',
    'NM'    : 	'Neisseria meningitidis',
    'SA'    : 	'Staphylococcus epidermidis',
    'SH'    : 	'Salmonella enterica',
    'SO'    : 	'Salmonella enterica',
    'SP'    : 	'Streptococcus pneumoniae',
    'ST'    : 	'Streptococcus pyogenes',
    'ST'    : 	'Streptococcus agalactiae',
    'V '    :   ' Vibrio cholerae',
    'CA'    :   ' Campylobacter'
}


create_dir(logs_dir)
logging.basicConfig(filename=logs_file, level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    force=True)

#############################################################
############################################################
# Start

get_info('[SD]','Script start running')

get_info('[SD]', 'Setup - preparing the data')
setup(genomes)

create_dir(specie_detection_dir)

arr_ecoli_shigella_reads, arr_others_reads = fetch_ecoli_and_shigella(genomes, ec_and_sh)

ec_sh_species   = gather_reads(arr_ecoli_shigella_reads)
others_species  = gather_reads(arr_others_reads)

get_info('[SD]', f'e.coli and shigella detected! -> {len(ec_sh_species)} specimen')
get_info('[SD]', f'other species -> {len(others_species)} specimen')


# Specie that are not ecoli shigella can be detected by kraken2 -- no problem!
if others_species:
    count = 0
    for sample, files in others_species.items():
        count += 1
        get_info('[SD]', f'Running kraken2 on {sample}      [{count}/{len(others_species)}]')
        r1 = files['R1']
        r2 = files['R2']
        
        kraken_report_dir = os.path.join(specie_detection_dir, 'kraken_report')
        create_dir(kraken_report_dir)
        report = f'{kraken_report_dir}/{sample}.report'
       
        metagenomics(
            r1=f'{genomes}/{r1}',
            r2=f'{genomes}/{r2}',
            report=report,
            db=DB
        )

        create_detection_report(
            input_file=report,
            sample=sample,
            output_file=detection_report_file
            )
        
get_info('[Report]','Report created!')
# If samples ECXXXX or F-ECXXXX or SGXXXXX are present, we need discriminate them
# in blasting genes specifi to ecoli and shigella, because kraken cannot
if ec_sh_species:
    count = 0
    get_info('[SD]', f'Detection of Ecoli, Shigella') 
   
    create_dir(ec_sg_dir)
    create_dir(assemblies_dir)

    for sample, files in ec_sh_species.items():
        get_info('[EC_SG]', f'Detection of {sample}      [{count}/{len(ec_sh_species)}]')

        r1 = files['R1']
        r2 = files['R2']

        trimmed_dir = os.path.join(samples_dir, sample, TRIMMED_DIR)
        create_dir(trimmed_dir)
        get_info('[EC_SG]', f'  -> Trimming')
        trimming(r1, r2, genomes, trimmed_dir)

        assembly_dir = os.path.join(samples_dir, sample, ASSEMBLY_DIR)
        create_dir(assembly_dir)
        get_info('[EC_SG]', f'  -> Assembly')
        assembly(
            r1=os.path.join(trimmed_dir, r1), 
            r2=os.path.join(trimmed_dir, r2),
            sample=sample,
            output=assembly_dir)
        
        if len(assembly_dir) > 0:
            gather_assemblies(
                source=f'{assembly_dir}/{sample}.fasta',
                destination=assemblies_dir
            )

    # Combine assemblies to one file to run blast on it
    combine_fastas(assemblies_dir, seq_combined_dir)
    
    create_dir(db_dir)
    create_dir(blast_dir)
    
    create_db(seq_combined_file, db_dir)

    for d in os.listdir(genes_dir):
        for f in os.listdir(f'{genes_dir}/{d}'):
            blast(
                query=os.path.join(genes_dir, d, f),
                db=f'{db_dir}/db',
                out=os.path.join(blast_dir, f'{f.split(".")[0]}.txt')
            )

    if os.path.exists(blast_dir):
        filter_blast_results(blast_dir,blast_report_file)
    
    # Check the specie - gene match
    if os.path.exists(blast_report_file):
        df = pd.read_csv(blast_report_file)
    
    values = []
    ec_founded = {}
    sg_founded = {}
    for index, row in df.iterrows():
        gene_name = row[GENE_NAME_COL]
        sample_name = row[SAMPLE_NAME_COL]

        # Since EIEC can have also ipaH, we need to check if the sample is in ec_founded
        # -> don't put it in the report as shigella too
        if sample_name in ec_founded:
            continue                
        
        value = species_by_genes[gene_name]

        if gene_name in ec_genes:
            ec_founded[sample_name] = True
        elif gene_name in sg_genes:
            sg_founded[sample_name] = True
        else:
            value = '?'
      
        values.append([sample_name, value])


    if values:
        header = [SAMPLE_NAME_COL, SPECIE_DETECTED_COL]
        df = pd.DataFrame(values, columns=header)
        df.to_csv(blast_sum_report_file, index=False)

    else:
        logging.error(f"[EC_SG Detection]   Error occurred while fetching report values\n")
        exit()



if os.path.exists(blast_sum_report_file):
    df_ec_sg = pd.read_csv(blast_sum_report_file)

    if os.path.exists(detection_report_file):
        get_info('[Report]', 'Combine reports')
        df_ec_sg.to_csv(detection_report_file, mode='a', index=False, header=False)
    else:
        df_ec_sg.to_csv(detection_report_file, index=False, header=False)


remove_dir(tmp_dir)   
remove_dir(samples_dir)  
# remove_dir(assemblies_dir)                    

# Create a separated file to get if matched or not
get_info('[Report]', 'Creating matched file')
df = pd.read_csv(detection_report_file)

no_matches = ''
for index, row in df.iterrows():
    sample           = row[SAMPLE_NAME_COL]
    specie_detected  = row[SPECIE_DETECTED_COL]
    sample_id        = re.split(r'(\d+)', sample)[0]
    specie_name      = species_list[sample_id].lower()
    
    if not species_list[sample_id].lower() == row[SPECIE_DETECTED_COL].lower():
        if len(no_matches) == 0:
            no_matches += "Species didn't matched:\n"
        no_matches += f"     -> {sample}: {specie_detected} detected instead of {specie_name}\n"

if len(no_matches) == 0:
    text = 'All the samples matched !'
else:
    text = no_matches

with open(matched_file, 'w') as file:
    file.write(text)







