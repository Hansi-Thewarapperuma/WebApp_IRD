'''
Function: Evaluation through genomic context criterion
Input: Dataframe
Output: A list of tuples containing index and evaluation (considerable_impact, no_considerable_impact)

Author: Hansi Thewarapperuma
date: 15/12/2023
'''

import pandas as pd

def evaluate_genomic_context(filtered_df):
    impact_column = 'ANN[0].IMPACT'
    ontology_terms_column = 'ANN[0].ANNOTATION'

    results = []

    for index, row in filtered_df.iterrows():
        # Check if the mentioned columns have non-null values
        if pd.notna(row[impact_column]) and pd.notna(row[ontology_terms_column]):
            if 'HIGH' in row[impact_column] or 'MODERATE' in row[impact_column]:
                # Check if 'ANN[0].ANNOTATION' is not one of the specified values
                if row[ontology_terms_column] not in ['synonymous_variant', 'intergenic_region', 'intron_variant',
                                                     'intragenic_variant', 'start_retained_variant', 'stop_retained_variant',
                                                     'initiater_codon_variant', '3_prime_UTR_variant', '5_prime_UTR_variant',
                                                     '5_prime_UTR_premature_start_codon_gain_variant',
                                                      ]:
                    results.append((index, 'considerable_impact'))
                else:
                    results.append((index, 'no_considerable_impact'))
            else:
                results.append((index, 'no_considerable_impact'))

    return results