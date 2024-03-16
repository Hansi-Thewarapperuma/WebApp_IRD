# web_app/app.py
import os
from flask import Flask, render_template, request, url_for
from MAF import evaluate_minor_allele_freq
from Clin_sig import evaluate_clinical_significance
from Cons_score import conservation_scores
from Func_pred import in_silico_functional_predictions
from Gen_context import evaluate_genomic_context
from file_process import tsv_to_df_filter_genes
from Quality import evaluate_quality

app = Flask(__name__)

# Set the upload folder
app.config['UPLOAD_FOLDER'] = 'uploads'
if not os.path.exists(app.config['UPLOAD_FOLDER']):
    os.makedirs(app.config['UPLOAD_FOLDER'])

# Limit the allowed extensions for file uploads
ALLOWED_EXTENSIONS = {'tsv'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# *********************** EVALUATION-RULES *************************************

# highest score 1 given to the highest priority evaluation value which ensures pathogenicity and vice versa
# All the evaluaton values were given the same priority

# 1 was assigned to the above terms as they support the pathogenicity of the variant
# -1 assigned to the opposite terms as they support the benignity of the variant
# 0 was assigned to the neutral terms which have no evidence about pathogenicity

# Define a scoring system globally      *********** VUS SCORE **********************
score_dict = {
    'considerable_impact': 1,
    'deleterious': 1,
    'rare': 1,
    'conserved': 1,
    'pathogenic': 1,
    'quality': 1,
    'vus': 0,


    'no_considerable_impact': -1,
    'not_deleterious': -1,
    'common': -1,
    'not_conserved': -1,
    'benign': -2,
    'less_quality': -1,

    'unresolved_clinical_significance': 0,
    'N/A': 0,
    'ambiguous_conservation': 0,
    'ambiguous_deleteriousness': 0
}


# Calculate the total score for each entry
def calculate_score(evaluations):
    total_score = sum(score_dict[eval_value] for eval_value in evaluations)
    return total_score


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process():

    try:
        # Check if the post request has the file part
        if 'user_tsv' not in request.files:
            return render_template('error.html', error='No file part')

        file = request.files['user_tsv']

        # If the user does not select a file, browser sends an empty file
        if file.filename == '':
            return render_template('error.html', error='No selected file')

        if file and allowed_file(file.filename):
            # Save the file to the upload folder
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)

            # Process the uploaded file (you can use this filepath in your tsv_to_df_filter_genes function)
            filtered_df = tsv_to_df_filter_genes(filepath)


    except Exception as e:
        return render_template('error.html', error=str(e))


    # user_tsv = request.form['user_tsv']
    # filtered_df = tsv_to_df_filter_genes(user_tsv)

    functional_pred = in_silico_functional_predictions(filtered_df)
    consv_score = conservation_scores(filtered_df)
    clinical_sig_scores = evaluate_clinical_significance(filtered_df)
    af_output = evaluate_minor_allele_freq(filtered_df)
    impact_output = evaluate_genomic_context(filtered_df)
    quality_status = evaluate_quality(filtered_df)

    results_dict = {}
    for result_list in [functional_pred, consv_score, clinical_sig_scores, af_output, impact_output, quality_status]:
        for result_tuple in result_list:
            index = result_tuple[0]
            evaluation = result_tuple[1:]

            if index in results_dict:
                results_dict[index].extend(evaluation)
            else:
                results_dict[index] = list(evaluation)


    # Filter out entries with 'N/A' as all the evaluation values
    # filtered_results_dict = {index: evaluations for index, evaluations in results_dict.items() if
    #                          'N/A' not in evaluations}

    # Sort entries based on total score
    sorted_results = sorted(results_dict.items(), key=lambda x: calculate_score(x[1]), reverse=True)

    # Print or use the sorted results as needed (********** COMMENT OUT TO PRINT ALL **************)
    # for index, evaluations in sorted_results:
    #     print(f"Index: {index}, Evaluations: {evaluations}, Total Score: {calculate_score(evaluations)}")

    # Find the highest total score
    if sorted_results:
        highest_total_score = calculate_score(sorted_results[0][1])
    else:
        highest_total_score = 0

    # Filter out entries with the highest total score
    filtered_highest_score_results = {index: evaluations for index, evaluations in sorted_results if
                                      calculate_score(evaluations) == highest_total_score}

    # Print or use the filtered results as needed : HIGHEST SCORE
    # for index, evaluations in filtered_highest_score_results.items():
    #     print(f"Index: {index}, Evaluations: {evaluations}, Total Score: {calculate_score(evaluations)}")

    #  ****** IN CASE IF YOU WANT TO CONSIDER THE SECOND HIGHEST SCORE*******
    # Find the second-highest total score
    if len(sorted_results) > 1:
        try:
            second_highest_total_score = calculate_score(next(
                evaluations for index, evaluations in sorted_results[1:] if
                calculate_score(evaluations) != highest_total_score))
        except StopIteration:
            second_highest_total_score = 0
    else:
        second_highest_total_score = 0

    # Filter out entries with the second-highest total score
    filtered_second_highest_score_results = {index: evaluations for index, evaluations in sorted_results if
                                             calculate_score(evaluations) == second_highest_total_score}

    # Print or use the filtered results as needed: SECOND HIGHEST SCORE
    # for index, evaluations in filtered_second_highest_score_results.items():
    #     print(f"Index: {index}, Evaluations: {evaluations}, Total Score: {calculate_score(evaluations)}")

    print(len(filtered_df))

    # Count the number of results with the highest score and second highest score
    pathogenic_count = len(filtered_highest_score_results)
    pathogenic_likely_pathogenic_count = len(filtered_second_highest_score_results)

    # Print the counts
    # print(f'Count of pathogenic variants detected: {pathogenic_count}')
    # print(f'Count of pathogenic or likely-pathogenic variants detected: {pathogenic_likely_pathogenic_count}')

    # Get the indexes
    pathogenic_indexes = list(filtered_highest_score_results.keys())
    pathogenic_likely_pathogenic_indexes = list(filtered_second_highest_score_results.keys())

    # Print the indexes with statements
    # if pathogenic_indexes:
    #     print('Consider these rows of your TSV for detected pathogenic variants:', pathogenic_indexes)

    # if pathogenic_likely_pathogenic_indexes:
    #     print('Consider these rows of your TSV for pathogenic or likely-pathogenic variants:',
    #           pathogenic_likely_pathogenic_indexes)

    # Write the detected pathogenic variants to a new TSV file
    def write_tsv_file(output_file, indexes):
        try:
            # Ensure indexes are within the valid range
            valid_indexes = [idx for idx in indexes if 0 <= idx < len(filtered_df)]

            if not valid_indexes:
                print("No valid indexes to write.")
                return

            output_df = filtered_df.iloc[valid_indexes]

            # Reset the index before saving to ensure consistent indexing in the saved TSV file
            output_df.reset_index(drop=True, inplace=True)

            output_df.to_csv(output_file, sep='\t', index=False, columns=filtered_df.columns)
            print(f'TSV file saved to {output_file}')
        except Exception as e:
            print(f'Error saving TSV file: {e}')


    # Specify the output file paths
    pathogenic_output_file = 'static/pathogenic_variants.tsv'
    pathogenic_likely_pathogenic_output_file = 'static/pathogenic_likely_pathogenic_variants.tsv'

    # Write the output TSV files
    if pathogenic_indexes:
        write_tsv_file(pathogenic_output_file, pathogenic_indexes)
        print(f'Pathogenic variants saved to {pathogenic_output_file}')

    if pathogenic_likely_pathogenic_indexes:
        write_tsv_file(pathogenic_likely_pathogenic_output_file, pathogenic_likely_pathogenic_indexes)
        print(f'Pathogenic or likely-pathogenic variants saved to {pathogenic_likely_pathogenic_output_file}')

        # Get data for pathogenic variants
        pathogenic_variants_data = []
        for index in pathogenic_indexes:
            variant_data = filtered_df.loc[index, ['CHROM', 'POS', 'ID', 'ANN[0].GENE']]
            pathogenic_variants_data.append(list(variant_data))

        # Get data for likely pathogenic variants
        likely_pathogenic_variants_data = []
        for index in pathogenic_likely_pathogenic_indexes:
            variant_data = filtered_df.loc[index, ['CHROM', 'POS', 'ID', 'ANN[0].GENE']]
            likely_pathogenic_variants_data.append(list(variant_data))

    # Render the results template with the data
    return render_template('results.html',
                           highest_score=filtered_highest_score_results,
                           second_highest_score=filtered_second_highest_score_results,
                           calculate_score=calculate_score,
                           pathogenic_count=pathogenic_count,
                           pathogenic_likely_pathogenic_count=pathogenic_likely_pathogenic_count,
                           pathogenic_indexes=pathogenic_indexes,
                           pathogenic_likely_pathogenic_indexes=pathogenic_likely_pathogenic_indexes,
                           pathogenic_output_file=url_for('static', filename='pathogenic_variants.tsv'),
                           pathogenic_likely_pathogenic_output_file=url_for('static', filename='pathogenic_likely_pathogenic_variants.tsv'),
                           pathogenic_variants_data=pathogenic_variants_data,
                           likely_pathogenic_variants_data=likely_pathogenic_variants_data

                           )



if __name__ == '__main__':
    app.run(debug=False)
