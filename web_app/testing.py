import pandas as pd

input_variants_df = pd.read_csv('CG408NEW_ann_extractedTSV.tsv', sep='\t', encoding='latin1')

# Read the virtual gene panel from the specified file
with open('web_app/gene_panel_file.txt', 'r') as gene_file:
    # Split each line into cancer type and associated genes
    virtual_gene_panel = {line.split(': ')[0]: line.split(': ')[1].split() for line in gene_file}

# Assuming gene names are in the column 'ANN[0].GENE', filter rows based on virtual gene panel
gene_column_name = 'ANN[0].GENE'
filtered_df = input_variants_df[
    input_variants_df[gene_column_name].isin([gene for genes in virtual_gene_panel.values() for gene in genes])]

# Resetting the index to retain the original index values
filtered_df.reset_index(drop=True, inplace=True)

print(filtered_df)