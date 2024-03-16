'''
Functon: TSV file conversion to dataframe
Input: TSV file
Output: Dataframe

Author: Hansi Thewarapperuma
Date: 25/12/2022
'''


# import pandas as pd
# def tsv_to_df_filter_genes(user_tsv):
#     # Define the virtual gene panel
#     virtual_gene_panel = {
#         "breast_ovarian_cancer": ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CHEK2", "EPCAM", "FANCC", "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "NF1", "PALB2", "PMS1", "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"],
#         "colorectal_cancer": ["APC", "ATM", "BMPR1A", "CDH1", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "MUTYH", "PMS2", "PTEN", "SMAD4", "STK11", "TP53"],
#         "thyroid_cancer": ["HRAS", "MUTYH", "PTEN", "RET", "SDHB", "SDHD", "TP53"],
#         "endometrial_cancer": ["BRCA1", "BRCA2", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "MUTYH", "PMS2", "PTEN",
#                                "TP53"],
#         "eye_cancer": ["BAP1", "RB1"],
#         "hereditary_diffuse_gastric_cancer": ["BMPR1A", "CDH1", "EPCAM", "KIT", "MLH1", "MSH2", "MSH6", "NF1", "PMS2",
#                                               "SDHB", "SDHC", "SDHD", "SMAD4", "STK11", "TP53"],
#         "adrenal_cancer": ["APC", "MEN1", "MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "TP53"],
#         "renal_urinary_tract_cancer": ["BAP1", "CDKN1C", "DICER1", "DIS3L2", "EPCAM", "FH", "FLCN", "GPC3", "MET",
#                                        "MLH1", "MSH2", "MSH6", "PMS1", "PTEN", "SDHB", "SDHC", "SMARCB1", "TP53",
#                                        "TSC1", "TSC2", "VHL", "WT1"],
#         "brain_neuro_system_cancer": ["AIP", "ALK", "APC", "DICER1", "EPCAM", "HRAS", "MEN1", "MLH1", "MSH2", "MSH6",
#                                       "NF1", "NF2", "PHOX2B", "PMS2", "PRKAR1A", "PTCH1", "PTEN", "RB1", "SMARCB1",
#                                       "SUFU", "TP53", "TSC1", "TSC2", "VHL"],
#         "liver_cancer": ["DICER1"],
#         "multiple_endocrine_neoplasia": ["MEN1"],
#         "familial_lung_cancer": ["BRCA1", "BRCA2", "CHEK2", "EGFR", "TP53"],
#         "pancreatic_cancer": ["APC", "ATM", "BMPR1A", "BRCA1", "BRCA2", "CDK4", "CDKN2A", "EPCAM", "FANCC", "MEN1",
#                               "MLH1", "MSH2", "MSH6", "NF1", "PALB2", "PMS2", "SMAD4", "STK11", "TP53", "TSC1", "TSC2",
#                               "VHL"]
#         # Add more cancer types and gene names as needed
#     }
#
#
#     try:
#         # Read the TSV file into a DataFrame using pandas
#         input_variants_df = pd.read_csv(user_tsv, sep='\t', encoding='latin1')
#
#         # Assuming gene names are in the column 'ANN[0].GENE', filter rows based on virtual gene panel
#         gene_column_name = 'ANN[0].GENE'
#         filtered_df = input_variants_df[input_variants_df[gene_column_name].isin([gene for genes in virtual_gene_panel.values() for gene in genes])]
#
#         # Resetting the index to retain the original index values
#         filtered_df.reset_index(drop=True, inplace=True)
#
#         return filtered_df
#     except Exception as e:
#         # Print an error message if reading the file fails
#         print(f"Error reading the file: {e}")
#         return None

import pandas as pd

def tsv_to_df_filter_genes(user_tsv):
    try:
        # Read the TSV file into a DataFrame using pandas
        input_variants_df = pd.read_csv(user_tsv, sep='\t', encoding='latin1')

        # Read the virtual gene panel from the specified file
        with open('gene_panel_file.txt', 'r') as gene_file:
            # Split each line into cancer type and associated genes
            virtual_gene_panel = {line.split(': ')[0]: line.split(': ')[1].split() for line in gene_file}

        # Assuming gene names are in the column 'ANN[0].GENE', filter rows based on virtual gene panel
        gene_column_name = 'ANN[0].GENE'
        filtered_df = input_variants_df[input_variants_df[gene_column_name].isin([gene for genes in virtual_gene_panel.values() for gene in genes])]

        # Resetting the index to retain the original index values
        filtered_df.reset_index(drop=True, inplace=True)

        return filtered_df
    except Exception as e:
        # Print an error message if reading the file fails
        print(f"Error reading the file: {e}")
        return None