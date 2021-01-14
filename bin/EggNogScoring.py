import pandas as pd
import numpy as np
import pickle
from Bio import SeqIO
import cobra
import warnings
import re

universal_model_path = 'data/reactionDatabase/bigModelv2.20b.sbml'
universal_csv_path = 'data/reactionDatabase/universal_v2.20.csv'

VERY_SMALL_NUMBER_NEG = -0.000000000001
VERY_SMALL_NUMBER = 0.00000000001
SCALE = 3
NEG_SCORE = -3


def get_gene_to_score_dict(fungi_id, annotation_file):
    """

    :param fungi_id: Fungi identification.
    :param annotation_file: eggNog file containing the functional annotations.
    :return: Dictionary containing the relation between ECs and the corresponding score, and dictionary containing
    the relation between the EC numbers and the corresponding genes.
    """
    # Opening the universal model information file to select the annotations of the EC numbers present in the
    # reaction database.
    universal = pd.read_csv(universal_csv_path)
    universal = universal.where((pd.notnull(universal)), "")

    # Collecting a list of all the EC numbers present in the universal
    ec_list = []
    for key, row in universal.iterrows():
        ECs = row['ECs'].split('|')
        for ec in ECs:
            ec_list.append(ec)
    ec_list = list(set(ec_list))
    ec_list = ec_list[1:len(ec_list)]

    # Opening the functional annotation file
    eggnog = pd.read_csv(annotation_file, skiprows=4, header=None, sep='\t')
    eggnog_selection = eggnog[eggnog[7].notnull()]

    genes = list(set(eggnog_selection[0].values))

    # getting the scores from the functional annotation file. It will be the maximum score of all the genes that
    # contains score for an EC number.

    # Dictionary containing dictionaries with the correspondence between the genes and the EC numbers.
    gene_to_ec_scores = {}

    for gene in genes:
        sel = eggnog[eggnog[0] == gene]
        ec_dictionary = {}
        for key, row in sel.iterrows():
            if gene in gene_to_ec_scores:
                ec_dictionary = gene_to_ec_scores[gene]
                ecs = row[7].split(',')
                for ec in ecs:
                    if ec in ec_dictionary:
                        ec_dictionary[ec] = np.max([ec_dictionary[row[7]], row[3]])
                    else:
                        ec_dictionary[ec] = row[3]
            else:
                ecs = row[7].split(',')
                for ec in ecs:
                    ec_dictionary[ec] = row[3]
            gene_to_ec_scores[gene] = ec_dictionary

    # goes through all the collected scores. If the EC number is partial it gets a very small negative number as score
    gene_to_ec_scores_corrected = {}
    for key, item in gene_to_ec_scores.items():
        ec_dictionary = {}
        for key2, item2 in item.items():
            new_item = item2
            if '.-' in key2[len(key2) - 2:len(key2)]:
                new_item = VERY_SMALL_NUMBER_NEG

            ec_dictionary[key2] = new_item
        gene_to_ec_scores_corrected[key] = ec_dictionary

    # Creates the dictionary EC to score and gets the maximum score for each EC number
    ec_to_score = {}
    for gene, ecs in gene_to_ec_scores_corrected.items():
        for ec, score in ecs.items():
            if ec in ec_to_score:
                score_list = ec_to_score[ec]
                score_list.append(score)
            else:
                ec_to_score[ec] = [score]

    ec_to_score_maximum = {}
    for ec, reaction_scores in ec_to_score.items():
        ec_to_score_maximum[ec] = np.max(reaction_scores)

    # Scaling the scores from 0 to SCALE
    avg_score = np.max(list(ec_to_score_maximum.values()))

    ec_to_score_final = {}
    for ec, score in ec_to_score_maximum.items():
        ec_to_score_final[ec] = score * SCALE / avg_score
        if score > avg_score:
            ec_to_score_final[ec] = SCALE

    # create EC to gene dictionary. Selecting only the genes that are near to the gene with maximum score.
    ecs_to_genes_dictionary = {}
    max_value = 0
    for ec in ec_list:
        genes = {}
        for key, item in gene_to_ec_scores_corrected.items():
            for key2, i in item.items():
                if ec == key2:
                    genes[key] = i
        if len(list(genes.values())) > 0:
            max_value = np.max(list(genes.values()))
        genes_list = []
        for key, item in genes.items():
            if item > max_value - 50:
                genes_list.append(str(key))
        ecs_to_genes_dictionary[ec] = list(set(genes_list))

    with open("results/" + fungi_id + "_GeneEc.npy", 'wb') as handle:
        pickle.dump(ecs_to_genes_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return ec_to_score_final, ecs_to_genes_dictionary


def assign_scores(ec_to_score_final):
    """
    Assings the scores collected to the ECs numbers to the corresponding reactions.
    :param ec_to_score_final: Dictionary with the scores corresponding to each EC number.
    :return: Pandas dataFrame of the universal csv containing a new column with the scores.
    """
    universal = pd.read_csv(universal_csv_path)
    universal = universal.where((pd.notnull(universal)), "")

    scores = []

    for key, row in universal.iterrows():
        universal_reaction_ec = row['ECs'].split('|')
        is_found = False
        sc = []
        for ec in universal_reaction_ec:
            ec_name = ec
            if len(ec) > 4 and ec_name in ec_to_score_final:
                sc.append(ec_to_score_final[ec_name])
                is_found = True

        if not is_found:
            if 'spontaneous' in row['ECs'] or 'artificial' in row['ECs']:
                scores.append(VERY_SMALL_NUMBER)
            elif 'transport' in row['ECs']:
                if '_e ' in row['Formula']:
                    scores.append(VERY_SMALL_NUMBER)
                else:
                    scores.append(VERY_SMALL_NUMBER)
            else:
                score_partial = get_partial_score(universal_reaction_ec)
                if score_partial == VERY_SMALL_NUMBER_NEG:
                    scores.append(score_partial)
                else:
                    scores.append(NEG_SCORE)

        else:
            scores.append(np.max(sc))

    universal['scores'] = scores

    return universal


def get_partial_score(universal_reaction_ec):
    """
    Set the score when the EC number of the reaction is partial.
    :param universal_reaction_ec: Reaction EC number.
    :return: Score.
    """
    for ec in universal_reaction_ec:
        if '.-' in ec:
            ec = ec[0:len(ec) - 2]
            return VERY_SMALL_NUMBER_NEG
        else:
            return NEG_SCORE


def get_compartment(metabolites):
    """
    From a given set of metabolites it returns the list of compartments where they belong.
    """
    comp = []
    for met in metabolites:
        comp.append(met.id[len(met.id) - 1])
    return comp


def get_genes(ec_numbers, ecs_genes_dictionary):
    """
    Function to get the genes associated to a EC number
    :param ec_numbers: String containing the list of EC numbers separated by "|".
    :param ecs_genes_dictionary: Dictionary of the ecs-genes association.
    :return: List of genes associated to a specific list of EC numbers.
    """
    try:
        genes = []
        ecs = ec_numbers.split('|')
        for ec in ecs:
            genes.append(ecs_genes_dictionary[ec])
        return genes
    except:
        return ""


def get_final_score(genes, compartment_label, cutoff, score, compartment_prediction_dataset):
    """
    Function to calculate the final reaction score by collecting all the compartment prediction of all the associated
    genes and using the maximum.
    :param genes: List of gene names to be checked.
    :param compartment_label: Label of the specific compartment prediction to be checked.
    :param cutoff: Cutoff used to decided the presence/absence of a reaction in a specific compartment.
    :param score: Reaction score obtained from the functional annotation.
    :param compartment_prediction_dataset: Pandas dataFrame containing the protein cellular localization prediction.
    :return: The score of a given reaction based on the compartment predictions.
    """
    compartment_prediction_score = []
    for gene in genes:
        for gen in gene:
            gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Id'] == gen]
            compartment_prediction_score.append(gen_selection.iloc[0][compartment_label])
    if len(compartment_prediction_score) > 0:
        compartment_prediction_score = np.max(compartment_prediction_score)
    else:
        return score
    return (0.00001 + (SCALE - score)) * (compartment_prediction_score - cutoff)


def scoring_compartments(sequence_file, compartment_file, universal, ecs_to_genes_dictionary, fungi_id):
    ids = []
    seq = []
    for record in SeqIO.parse(sequence_file, "fasta"):
        ids.append(str(record.id))
        seq.append(str(record.seq))

    compartment_predictions = np.load(compartment_file)

    endo = []
    mito = []
    pero = []
    other = []

    for pred in compartment_predictions:
        endo.append(pred[0])
        mito.append(pred[1])
        pero.append(pred[2])
        other.append(pred[3])

    dataset = pd.DataFrame(
        {'Id': ids, 'Sequences': seq, 'Endoplasmic reticulum': endo, 'Mitochondria': mito, 'Peroxisome': pero,
         'Other': other})

    # calculate final scores for each reaction in the universal model
    model = cobra.io.read_sbml_model(universal_model_path)

    final_scores = {}
    for reaction in model.reactions:
        if reaction.id.endswith('_E'):
            continue
        id = re.search('UF\d{5}', reaction.id)
        try:
            id = id.group(0)
        except:
            id = "BIOMASS"

        try:
            sel = universal[universal['IDs'].str.contains(id)]
            if len(sel) > 0:
                final_scores['R_' + reaction.id] = sel.iloc[0]['scores']
            if '_T' in reaction.id or '_t' in reaction.id:
                final_scores['R_' + reaction.id] = VERY_SMALL_NUMBER
        except:
            final_scores['R_' + reaction.id] = VERY_SMALL_NUMBER_NEG

    final_scores_with_compartments = final_scores.copy()
    for reaction in model.reactions:
        name = ""
        if reaction.id.startswith('R_'):
            name = reaction.id
        else:
            name = 'R_' + reaction.id
        if reaction.id.endswith('_E'):
            continue

        comp = get_compartment(reaction.metabolites)
        comp = list(set(comp))
        # print(reaction.id,comp)
        if len(comp) == 1:
            if reaction.id.startswith('EX'):
                continue
            name2 = re.search('UF\d{5}', reaction.id)
            try:
                name2 = name2.group(0)
            except:
                name2 = "BIOMASS"

            sel = universal[universal['IDs'].str.contains(name2)]
            if len(sel) > 0:
                score = final_scores[name]
                if score < 0:
                    final_scores_with_compartments[name] = score
                else:
                    ec_number = sel.iloc[0]['ECs']
                    genes = get_genes(ec_number, ecs_to_genes_dictionary)
                    compartment_score=0
                    if 'm' in comp:
                        compartment_score = get_final_score(genes, "Mitochondria", 0.2, score, dataset)
                    if 'x' in comp:
                        compartment_score = get_final_score(genes, "Peroxisome", 0.15, score, dataset)
                    if 'c' in comp:
                        compartment_score = get_final_score(genes, "Other", 0.2, score, dataset)
                    if 'r' in comp:
                        compartment_score = get_final_score(genes, "Endoplasmic reticulum", 0.15, score, dataset)
                    if 'n' in comp:
                        compartment_score = get_final_score(genes, "Other", 0.2, score, dataset)
                    if 'g' in comp:
                        compartment_score = get_final_score(genes, "Other", 0.15, score, dataset)

                    if compartment_score < 0:
                        final_scores_with_compartments[name] = compartment_score
                    else:
                        final_scores_with_compartments[name] = score

    real_ids = []
    ids = final_scores.keys()
    for id in ids:
        name = re.search('UF\d{5}', id)
        try:
            name = name.group(0)
        except:
            name = "BIOMASS"
        real_ids.append(name)
    real_ids = list(set(real_ids))

    def getIds(id):
        ids = []
        for key, score in final_scores.items():
            if id in key:
                ids.append(key)
        return ids

    for id in real_ids:
        ids = getIds(id)
        anyPos = False
        for rxnId in ids:
            if final_scores[rxnId] > 0:
                anyPos = True
                break
        if anyPos:
            allNeg = True
            for rxnId in ids:
                if final_scores_with_compartments[rxnId] > 0:
                    allNeg = False
                    break
            if allNeg:
                for rxnId in ids:
                    final_scores_with_compartments[rxnId] = final_scores[rxnId]

    # Adding artificial scores to these reactions (CO2, water and O2 exchanges and complex IV"
    final_scores_with_compartments["R_UF03227_TCE"] = 6
    final_scores_with_compartments["R_UF00806_M"] = 6
    final_scores_with_compartments["R_UF03314_TCE"] = 6
    final_scores_with_compartments["R_UF03382_TCE"] = 6


    # Creating a dataframe with the final scores for reporting
    ids = []
    score = []
    formula = []
    definition = []
    name = []
    ec = []
    path = []
    annotation = pd.DataFrame()
    for key, val in final_scores_with_compartments.items():
        id = re.search('UF\d{5}', key)
        try:
            id = id.group(0)
            sel = universal[universal['IDs'].str.contains(id)]
            reaction = model.reactions.get_by_id(key[2:len(key)])
            if len(sel) > 0:
                formula.append(reaction.build_reaction_string())
                definition.append(sel.iloc[0]['Definition'])
                ec.append(sel.iloc[0]['ECs'])
                name.append(sel.iloc[0]['Reaction name'])
                path.append(sel.iloc[0]['Pathways'])
                score.append(val)
                ids.append(key)
            else:
                print(key)
        except:
            reaction = model.reactions.get_by_id(key[2:len(key)])
            formula.append(reaction.build_reaction_string())
            definition.append(reaction.build_reaction_string(use_metabolite_names=True))
            ec.append("")
            name.append(reaction.name)
            path.append("")
            score.append(val)
            ids.append(key)
            # print(key)
    annotation['IDs'] = ids
    annotation['Score'] = score
    annotation['Reaction Name'] = name
    annotation['Formula'] = formula
    annotation['Definition'] = definition
    annotation['EC'] = ec
    annotation['Pathways'] = path
    annotation.columns = ['IDs', 'Score', 'Reaction Name', 'Formula', 'Definition', 'EC', 'Pathways']
    annotation.to_csv("results/" + fungi_id + "_annotations.csv", index=False)

    print("Scores calculated")

    return final_scores_with_compartments
