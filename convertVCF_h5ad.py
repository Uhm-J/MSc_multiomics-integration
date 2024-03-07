import argparse
import gzip
import pandas as pd
import anndata as ad
import scipy.sparse as sp


def parse_genotype(genotype_str):
    """Parse genotype string to a numerical format."""
    genotype_str = genotype_str.split(':')[0]
    if pd.isna(genotype_str):
        return  0  # Set missing values to 0
    if genotype_str == '0/1':
        return 2  # Heterozygous
    if genotype_str == '1/1':
        return 3  # Homozygous variant
    return 1  # Homozygous reference or any other case

def get_annotation(vcf_data):
    """Preprocess VCF data to extract genotype and phenotype information."""
    print("Preprocessing Annotation data...")
    n_snvs = len(vcf_data)

    # Function to extract AF, AR2, and DR2 values
    def extract_info_values(info):
        info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
        af = float(info_dict.get('AF', 0))
        ar2 = float(info_dict.get('AR2', 0))
        dr2 = float(info_dict.get('DR2', 0))
        return af, ar2, dr2
    
    
    # Apply the extraction function to each row in the INFO column
    extracted_values = vcf_data['INFO'].apply(extract_info_values)
    
    # Split the tuple into separate columns
    df = pd.DataFrame(extracted_values.tolist(), columns=['AF', 'AR2', 'DR2'])

    df['ref_allele'] = vcf_data['REF']
    df['alt_allele'] = vcf_data['ALT']
    df['chrom'] = vcf_data['CHROM']
    df['pos'] = vcf_data['POS']
    df['id'] = df['chrom'] + ':' + df['pos'].astype(str) + '_' + df['ref_allele'] + '/' + df['alt_allele']

    print(df)

    print("Preprocessing complete.")
    
    return df

def get_cell_annotation(vcf_data):
    cell_names = vcf_data.columns[9:]
    phenotype_labels = pd.Series([
        "Bipolar" if "Bipolar" in name 
        else "Tripolar" if "Tripolar" in name 
        else "Control" if 'Control' in name 
        else 'NaN' for name in cell_names
    ])
    batch_labels = pd.Series([c.split('_')[0] for c in cell_names]) if args.batch is None else pd.Series([args.batch for c in cell_names])
    return pd.DataFrame({'Phenotype': phenotype_labels, 'Batch': batch_labels})


def get_genotype(vcf_data):
    """Preprocess VCF data to extract genotype and phenotype information."""
    print("Preprocessing VCF data...")
    n_snvs = len(vcf_data)

    cell_names = vcf_data.columns[9:]
    
    print(vcf_data[cell_names])
    genotype_data = vcf_data[cell_names]
    genotype_matrix = genotype_data.map(parse_genotype)
    print(genotype_matrix)
    print("Preprocessing complete.")
    
    # Transpose the matrix
    transposed_genotype_matrix = genotype_matrix.T

    # Convert the transposed matrix to CSR format
    csr_matrix = sp.csr_matrix(transposed_genotype_matrix)

    return csr_matrix, cell_names

def read_vcf(vcf_path):
    """Read VCF file, correctly handling metadata lines and the header."""
    print(f"Reading VCF file: {vcf_path}")
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    with open_func(vcf_path, 'rt') as file:
        for line in file:
            if line.startswith('#CHROM'):
                header_line = line.strip()
                break
    columns = header_line.lstrip('#').split('\t')
    vcf_data = pd.read_csv(vcf_path, comment='#', sep='\t', names=columns, low_memory=False)
    return vcf_data

def main():
    parser = argparse.ArgumentParser(description="TBD")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the VCF file.")
    parser.add_argument("--out", type=str, default="germline_snvs", help="Path to the VCF file.")
    parser.add_argument("--batch", type=str, help="Path to the VCF file.")
    global args

    args = parser.parse_args()
    print(f"Opening {args.vcf}")
    vcf_data = read_vcf(args.vcf)
    print(vcf_data.columns.tolist)
    sparse_genotype_matrix, cell_ids = get_genotype(vcf_data)
    snv_annotation = get_annotation(vcf_data)
    cell_annotation = get_cell_annotation(vcf_data)


    snv_ids = snv_annotation.pop('id')
    
    if args.batch is not None:
        cell_ids = [args.batch + '_' + c for c in cell_ids]
    print(cell_annotation) 
    print(type(cell_annotation))

    # Create AnnData object from genotype matrix and annotations
    adata = ad.AnnData(sparse_genotype_matrix)  # transpose genotype matrix to match anndata convention
    adata.var = snv_annotation
    adata.var_names = snv_ids
    adata.obs = cell_annotation
    adata.obs_names = cell_ids


    print(adata.X)
    print(adata.var)
    print(adata.obs)
    print(adata.shape)

    adata.write_h5ad(f"{args.out}.h5ad")

if __name__ == "__main__":
    main()