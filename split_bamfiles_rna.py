import pysam
import os
import argparse
import re
from tqdm import tqdm
import subprocess
import pandas as pd

def add_or_replace_read_groups(input_bam, output_bam, cell_id, barcode, sample, phenotype=''):
    # Construct the Picard command
    output_bam = output_bam.split('/')
    output_dir = '/'.join(output_bam[:-1])
    output_name = output_bam[-1].split('.')[0] + f"_ReadGroup.bam"
    output_file = f"{output_dir}/{output_name}"
    picard_command = [
        "java", "-jar", "/home/jorrit/picard.jar", "AddOrReplaceReadGroups",
        "I=" + input_bam,
        "O=" + output_file,
        "RGID="+ sample + "_" + cell_id,
        "RGLB=SORT-seq",
        "RGPL=SCD",
        "RGPU=" + barcode,
        "RGSM=" + cell_id
    ]

    # Run the Picard command
    subprocess.run(picard_command)

def read_phenotypes_from_excel(file_path: str, skip_rows: int) -> pd.DataFrame:
    # Read the specified range
    data = pd.read_excel(file_path, engine='openpyxl', usecols="A:Y", skiprows=range(0, skip_rows), nrows=17)
    pt = pd.read_excel(file_path, engine='openpyxl', usecols="AC:AE", skiprows=range(0, skip_rows+1), nrows=8)
    pt.columns = ['value', 'phenotype', 'phenotype2']
    pt = pt.set_index('value')['phenotype'].to_dict()
    data.replace(pt, inplace=True)
    #data = data.set_index('Unnamed: 0')
    # Use melt to transform the DataFrame
    # Assuming 'Unnamed: 0' is the column with row letters (A, B, C, etc.)
    melted_data = pd.melt(data, id_vars=['Unnamed: 0'], var_name='Column', value_name='Phenotype')
    
    # Creating a combined identifier as you wanted, e.g., A1, A2, etc.
    melted_data['Identifier'] = melted_data['Unnamed: 0'] + melted_data['Column'].astype(str)
    
    # Selecting only the combined identifier and phenotype columns
    final_data = melted_data[['Identifier', 'Phenotype']].dropna()
    
    return final_data

def create_well_annotations(phenotype_mapping, well_barcode_match_file):
    """Create a dictionary with well annotations from the well barcode match file."""
    well_annotations = {}
    with open(well_barcode_match_file, 'r') as file:
        for line in file:
            well,_, barcode = line.strip().split('\t')
            well_annotations[barcode] = well

    barcode_df = pd.DataFrame(list(well_annotations.items()), columns=['Barcode', 'Identifier'])
    print(barcode_df)
    combined_df = pd.merge(phenotype_mapping, barcode_df, on='Identifier')

    return combined_df

def extract_info(bc, annotation):
    result = annotation[annotation['Barcode'] == bc]
    
    # Check if there are any matches
    if not result.empty:
        # Assuming each barcode is unique and thus only one match is found
        info = result.iloc[0]  # Get the first (and presumably only) match
        return {
            'Identifier': info['Identifier'],
            'Phenotype': info['Phenotype']
        }
    else:
        return None


def split_bam_by_cell(input_bam, output_dir, sample):
    # Open the input BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")

    # Dictionary to store open file handles for each cell
    cell_files = {}
    unknown_barcodes = set()
    # Iterate over all reads in the BAM file using tqdm
    for read in tqdm(bamfile, desc="Splitting BAM file"):
        # Extract the cell barcode (CR tag)
        try:
            cell_barcode = read.get_tag('CB')
            cell = extract_info(cell_barcode, well_annotations)

        except KeyError:
            # Skip reads without a CR tag
            continue
        except:
            print("Unexpected error:", sys.exc_info()[0])
            continue
        
        # Use the cell barcode directly as the unique identifier
        cell = extract_info(cell_barcode, well_annotations)
        if cell_barcode in unknown_barcodes:
            continue
        
        if cell is None:
            print(f"Skipping read with barcode {cell_barcode}")
            unknown_barcodes.add(cell_barcode)
            continue
        
        unique_id = f"{cell['Identifier']}_{cell_barcode}"
        # Open a new file for this cell if it doesn't exist yet
        if unique_id not in cell_files:
            print(f"Creating file for cell {unique_id} with name {sample}_bi_{cell['Identifier']}_{cell['Phenotype']}.bam")
            cell_filename = os.path.join(output_dir, f"{sample}_bi_{cell['Identifier']}_{cell['Phenotype']}.bam")
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
    print(f"Unknown barcodes: {unknown_barcodes}")
    print(f"Number of unknown barcodes: {len(unknown_barcodes)}")

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
    
    # Read phenotypes from the annotation file

    phenotype_mapping = read_phenotypes_from_excel(args.annotation_file, args.skip_rows)
    print(phenotype_mapping)
    well_barcode_match_file = "/mnt/d/JorritvU/Tripolar/scRNA-seq/SORTseq_cellbarcodes.tsv"
    global well_annotations

    well_annotations = create_well_annotations(phenotype_mapping, well_barcode_match_file)
    print(well_annotations)
    print(extract_info("CGTGATTC", well_annotations))

    # Split the BAM file by cell
    split_bam_by_cell(args.input, args.output_dir, args.sample)

    # Rename the files in the output directory

if __name__ == "__main__":
    main()