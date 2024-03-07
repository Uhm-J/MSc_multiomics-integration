import pysam
import os
import argparse
import re
from tqdm import tqdm
import subprocess
import pandas as pd

def extract_tags(read_name):
    # Regular expressions to match the cell identifier (bi tag) and barcode (bc tag)
    bi_match = re.search(r'bi:(\d+)', read_name)
    bc_match = re.search(r'bc:([A-Za-z0-9]+)', read_name)

    cell_id = bi_match.group(1) if bi_match else None
    barcode = bc_match.group(1) if bc_match else None

    return cell_id, barcode

def add_or_replace_read_groups(input_bam, output_bam, cell_id, barcode, sample):
    # Construct the Picard command
    output_bam = output_bam.split('/')
    output_dir = '/'.join(output_bam[:-1])
    output_name = f"{sample}_bi_{cell_id}.bam"
    output_file = f"{output_dir}/{output_name}"
    picard_command = [
        "java", "-jar", "/home/jorrit/picard.jar", "AddOrReplaceReadGroups",
        "I=" + input_bam,
        "O=" + output_file,
        "RGID="+ sample + "_" + cell_id,
        "RGLB=NlaIII",
        "RGPL=SCC",
        "RGPU=" + barcode,
        "RGSM=" + cell_id
    ]

    # Run the Picard command
    subprocess.run(picard_command)

def split_bam_by_cell(input_bam, output_dir, sample):
    # Open the input BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")

    # Dictionary to store open file handles for each cell
    cell_files = {}

    # Iterate over all reads in the BAM file using tqdm
    for read in tqdm(bamfile, desc="Splitting BAM file"):
        # Extract the cell identifier from the read name
        cell_id = None       
        if read.query_name:
            cell_id, barcode = extract_tags(read.query_name)
        # Skip reads without a cell identifier
        if cell_id is None:
            continue

        unique_id = f"{cell_id}_{barcode}"
        # Open a new file for this cell if it doesn't exist yet
        if unique_id not in cell_files:
            print(f"Creating file for cell {unique_id}")
            cell_filename = os.path.join(output_dir, f"{sample}_bi_{cell_id}-temp-.bam")
            cell_files[unique_id] = pysam.AlignmentFile(cell_filename, "wb", template=bamfile)

        # Write read to the appropriate file
        cell_files[unique_id].write(read)

    # Close all open files
    for unique_id, file in cell_files.items():
        file.close()
        cell_id, barcode = unique_id.split("_")
        # Add read groups to the BAM file
        add_or_replace_read_groups(file.filename.decode(), file.filename.decode(), cell_id, barcode, sample)
    bamfile.close()

def read_phenotypes_from_excel(file_path: str, skip_rows: int) -> dict:
    # Read the specified range
    data = pd.read_excel(file_path, engine='openpyxl', usecols="A:Y", skiprows=range(0, skip_rows), nrows=17)
    pt = pd.read_excel(file_path, engine='openpyxl', usecols="AC:AE", skiprows=range(0, skip_rows+1), nrows=8)
    pt.columns = ['value', 'phenotype', 'phenotype2']
    pt = pt.set_index('value')['phenotype'].to_dict()
    pt[0] = "Empty"
    data.replace(pt, inplace=True)
    data = data.set_index('Unnamed: 0')
    
    # Flatten the dataframe to create a mapping of cell_id to phenotype
    phenotype_mapping = {}
    cell_id = 0
    for index, row in data.iterrows():
        for col, value in row.items():
            cell_id += 1
            identifier = f"bi_{cell_id}"
            phenotype_mapping[identifier] = value
    return phenotype_mapping

def rename_files(phenotype_mapping: dict, directory: str, sample: str):
    # Now iterate over the files in the directory
    for filename in os.listdir(directory):
        # Split the filename into the sample, cellID, and extension
        if "-temp-" in filename:
            continue
        parts = filename.split(".")
        extension = '.'.join(parts[1:])
        parts = parts[0].split('_')
        sample = parts[0]
        cell_id = '_'.join(parts[1:]).split('-')[0]  # Join back all parts except the last one
        
        # If the cell_id is in our mapping, rename the file
        if cell_id in phenotype_mapping:
            new_filename = f"{sample}_{cell_id}_{phenotype_mapping[cell_id]}.{extension}"
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
        else:
            print(f"No phenotype mapping for {cell_id}, file {filename} not renamed.")



def main():
    parser = argparse.ArgumentParser(description='Split BAM file by cell identifier.')
    parser.add_argument('--input', required=True, help='Input BAM file')
    parser.add_argument('--output_dir', required=True, help='Output directory for split BAM files')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--annotation_file', required=True, help='Excel file with annotations')
    parser.add_argument('--skip_rows', type=int, required=True, help='Number of rows to skip in the annotation file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    split_bam_by_cell(args.input, args.output_dir, args.sample)

    # Read phenotypes from the annotation file
    phenotype_mapping = read_phenotypes_from_excel(args.annotation_file, args.skip_rows)
    print(phenotype_mapping)
    # Rename the files in the output directory
    rename_files(phenotype_mapping, args.output_dir, args.sample)

if __name__ == "__main__":
    main()