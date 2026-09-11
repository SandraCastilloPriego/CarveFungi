import pandas as pd
import numpy as np
import pickle
import cobra
import re
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(script_dir, '..', 'data')
data_path = os.path.normpath(data_path)
results_path = os.path.join(script_dir, '..', 'results')
results_path = os.path.normpath(results_path)

universal_model_path = os.path.join(data_path, 'reactionDatabase/bigModelv3.8b.sbml')
universal_csv_path = os.path.join(data_path, 'reactionDatabase/universal_v3.8.csv')

VERY_SMALL_NUMBER_NEG = -0.000000000001
VERY_SMALL_NUMBER = 0.00000000001
SCALE = 3
NEG_SCORE = -3


def _read_annotation_with_dynamic_header(path):
    """Read a CSV/TSV file, auto-detecting delimiter and header location."""
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if '#query' in line:
                # Detect delimiter by checking the header line
                if '\t' in line:
                    sep = '\t'
                elif ',' in line:
                    sep = ','
                else:
                    sep = '\t'  # default to tab if unclear
                return pd.read_csv(path, skiprows=i, header=0, sep=sep)
    raise ValueError("Could not find '#query' column in file")



def _extract_universal_ecs():
    """Return a set of all EC numbers appearing in the universal reaction database."""
    uni = pd.read_csv(universal_csv_path).fillna('')
    ec_sets = (row.split('|') for row in uni['ECs'])
    ecs = {ec for group in ec_sets for ec in group}
    # original implementation dropped the first element for some reason; keep behaviour
    ecs.discard('')
    return ecs


def _score_by_gene(annotations):
    """Build a nested dict gene->ec->max_score from the annotation dataframe."""
    gene_ec = {}
    for _, row in annotations.iterrows():
        gene = row['#query']
        score = row['score']
        for ec in row['EC'].split(','):
            gene_ec.setdefault(gene, {}).setdefault(ec, score)
            gene_ec[gene][ec] = max(gene_ec[gene][ec], score)
    return gene_ec


def _correct_partial_ecs(gene_ec):
    """Replace scores for partial ECs (contains '.-') with a tiny negative value."""
    out = {}
    for gene, ecs in gene_ec.items():
        out[gene] = {}
        for ec, sc in ecs.items():
            if ec.endswith('.-'):
                out[gene][ec] = VERY_SMALL_NUMBER_NEG
            else:
                out[gene][ec] = sc
    return out


def _compute_ec_max_scores(gene_ec):
    """From gene->ec->score compute ec->maximum score across all genes."""
    ec_scores = {}
    for ecs in gene_ec.values():
        for ec, sc in ecs.items():
            ec_scores.setdefault(ec, []).append(sc)
    return {ec: np.max(vals) for ec, vals in ec_scores.items()}


def _scale_scores(ec_max):
    """Scale scores so that the highest value becomes SCALE (linear)."""
    if not ec_max:
        return {}
    max_val = np.max(list(ec_max.values()))
    scaled = {}
    for ec, sc in ec_max.items():
        val = sc * SCALE / max_val
        scaled[ec] = SCALE if sc > max_val else val
    return scaled


def _select_genes_for_ec(ec_list, gene_ec, compartment_dataset):
    """Produce the EC->[[gene, compartments],...] mapping used downstream."""
    mapping = {}
    for ec in ec_list:
        scores = {g: ecs[ec] for g, ecs in gene_ec.items() if ec in ecs}
        if not scores:
            mapping[ec] = []
            continue
        max_score = np.max(list(scores.values()))
        candidates = [g for g, sc in scores.items() if sc > max_score - 120]
        genes_with_comp = []
        for gene in set(candidates):
            comp_list = []
            
            comp_series = _get_compartment_series(compartment_dataset, gene)
            #print(gene, comp_series)
            if isinstance(comp_series, pd.Series):
                for label, thresh in [('M', 0.31), ('P', 0.04), ('O', 0.40),
                                       ('E', 0.29), ('N', 0.23)]:
                    if comp_series.get(label, 0) > thresh:
                        comp_list.append(label)
            if not comp_list:
                comp_list = ['C']
            genes_with_comp.append([gene, comp_list])
        mapping[ec] = genes_with_comp
    return mapping


def get_gene_to_score_dict(fungi_id, annotation_file, compartment_prediction_dataset, results_path=None):
    """Compute EC scores and associated genes from eggNOG annotations.

    **Returns**: ``(ec_to_score_final, ecs_to_genes_dictionary)``

    - ``ec_to_score_final`` is an EC->scaled score dictionary (0..SCALE).
    - ``ecs_to_genes_dictionary`` is an EC->list of ``[gene, [compartments]]`` pairs.

    The computation performs the following steps:
    1. load universal model and gather all EC numbers
    2. read the annotation table and collapse scores by gene/EC
    3. correct partial ECs, take maximum score per EC and rescale
    4. select genes near the top score for each EC and determine
       their compartment labels using the provided dataset
    5. persist the gene/EC mapping in ``results_path``
    """
    if results_path is None:
        results_path = globals()['results_path']
    
    # Ensure results directory exists
    os.makedirs(results_path, exist_ok=True)
    
    universal_ecs = _extract_universal_ecs()
    annotations = _read_annotation_with_dynamic_header(annotation_file)
    annotations = annotations[annotations['EC'].notnull()]
    

    gene_ec = _score_by_gene(annotations)
    gene_ec = _correct_partial_ecs(gene_ec)

    ec_max = _compute_ec_max_scores(gene_ec)
    ec_to_score_final = _scale_scores(ec_max)

    ecs_to_genes_dictionary = _select_genes_for_ec(
        universal_ecs, gene_ec, compartment_prediction_dataset)

    gene_ec_path = os.path.join(results_path, f"{fungi_id}_GeneEc.npy")
    with open(gene_ec_path, 'wb') as handle:
        pickle.dump(ecs_to_genes_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
    #print(f"DEBUG: Saved GeneEc file to {gene_ec_path}")

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


def _get_compartment_series(compartment_prediction_dataset, gene):
    """Return a pandas Series-like mapping for a given gene from the
    compartment prediction dataset. The dataset may be a DataFrame, a
    dict mapping gene->row, or a path to a CSV file. Returns None when the
    gene is not found.
    """
    # If a path was passed, try to read it
    if isinstance(compartment_prediction_dataset, str):
        try:
            #print(f"DEBUG: Reading compartment prediction file: {compartment_prediction_dataset}")
            # Try with default comma separator first
            compartment_prediction_dataset = pd.read_csv(compartment_prediction_dataset)
            #print(f"DEBUG: Successfully read compartment prediction file: {len(compartment_prediction_dataset)} rows")
        except Exception as e:
            #print(f"DEBUG: Failed to read compartment prediction file with comma separator: {e}")
            try:
                # Try with tab separator
                compartment_prediction_dataset = pd.read_csv(compartment_prediction_dataset, sep='\t')
                #print(f"DEBUG: Successfully read with tab separator")
            except Exception as e2:
                #print(f"DEBUG: Failed to read with tab separator: {e2}")
                return None

    # DataFrame case
    if isinstance(compartment_prediction_dataset, pd.DataFrame):
        try:
            if "." in gene:
                gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'] == gene]
            else:               
                gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'] == gene]
                if len(gen_selection) == 0:
                    gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'].str.contains(gene)]
           
            if len(gen_selection) > 0:
                return gen_selection.iloc[0]
            return None
        except Exception:
            return None

    # dict-like case
    if isinstance(compartment_prediction_dataset, dict):
        val = compartment_prediction_dataset.get(gene)
        if val is None:
            return None
        # Convert dict/list to Series for uniform access
        try:
            return pd.Series(val)
        except Exception:
            return None

    return None


def get_compartment(metabolites):
    """
    From a given set of metabolites it returns the list of compartments where they belong.
    """
    comp = []
    for met in metabolites:
        comp.append(met.id[len(met.id) - 1])
    return comp


def get_genes(ec_numbers, ecs_genes_dictionary):
    genes = []
    ecs = ec_numbers.split('|')
    for ec in ecs:
        if ec in ecs_genes_dictionary:  # gracefully skips missing ECs
            genes.append(ecs_genes_dictionary[ec])
    return genes if genes else ""

def flatten_list(mixed_list):
    # If it's already a flat list, return it as-is
    if all(isinstance(item, str) for item in mixed_list):
        return mixed_list
    # Otherwise, flatten the nested list
    return [item for sublist in mixed_list for item in sublist]


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
    gene = flatten_list(genes)
    #print("Genes:", gene)
    for gen in gene:
        #print("Gene:", gen)
        gen = gen[0]
        try:
            #print(gen, compartment_label, gen.split(".")[0])
            if "." in gen:
                gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'] == gen]
            else:               
                gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'] == gen]
                if len(gen_selection) == 0:
                    gen_selection = compartment_prediction_dataset[compartment_prediction_dataset['Ids'].str.contains(gen)]
            compartment_prediction_score.append(gen_selection.iloc[0][compartment_label])
        except:
            pass
    
    if len(compartment_prediction_score) > 0:
        compartment_prediction_score = np.max(compartment_prediction_score)
    else:
        return score
    '''if compartment_prediction_score < cutoff:
        return NEG_SCORE
    else:
        return score'''
    return (1 + (SCALE - score)) * (compartment_prediction_score - cutoff)


def scoring_compartments(compartment_file, universal, ecs_to_genes_dictionary, fungi_id, results_path=None):

    if results_path is None:
        results_path = globals()['results_path']
    
    # Ensure results directory exists
    os.makedirs(results_path, exist_ok=True)

    # Read compartment file with better error handling
    #print(f"DEBUG: Reading compartment file: {compartment_file}")
    try:
        # Try with default comma separator first
        dataset = pd.read_csv(compartment_file)
        #print(f"DEBUG: Successfully read compartment file with comma separator: {len(dataset)} rows")
    except Exception as e:
        #print(f"DEBUG: Failed with comma separator: {e}")
        try:
            # Try with tab separator
            dataset = pd.read_csv(compartment_file, sep='\t')
            #print(f"DEBUG: Successfully read compartment file with tab separator: {len(dataset)} rows")
        except Exception as e2:
            #print(f"DEBUG: Failed with tab separator: {e2}")
            raise ValueError(f"Could not read compartment file {compartment_file}: {e}")

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
                final_scores[reaction.id] = sel.iloc[0]['scores']
            else:
                final_scores[reaction.id] = VERY_SMALL_NUMBER
            if '_T' in reaction.id or '_t' in reaction.id:
                final_scores[reaction.id] = VERY_SMALL_NUMBER
        except:
            final_scores[reaction.id] = VERY_SMALL_NUMBER_NEG

    final_scores_with_compartments = final_scores.copy()
    for reaction in model.reactions:          
        name = reaction.id
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
                if score < 0: #or sel.iloc[0]['Compartments'] == 0:
                    final_scores_with_compartments[name] = score
                else:
                    ec_number = sel.iloc[0]['ECs']
                    genes = get_genes(ec_number, ecs_to_genes_dictionary)
                    print(name, genes, comp)
                    compartment_score=0
                    if 'm' in comp:
                        compartment_score = get_final_score(genes, "M", 0.31, score, dataset)
                    if 'x' in comp:
                        compartment_score = get_final_score(genes, "P", 0.04, score, dataset)
                    if 'c' in comp:
                        compartment_score = get_final_score(genes, "O", 0.40, score, dataset)
                    if 'r' in comp:
                        compartment_score = get_final_score(genes, "E", 0.29, score, dataset)
                    if 'n' in comp:
                        compartment_score = get_final_score(genes, "N", 0.23, score, dataset)
                    if 'l' in comp:
                            compartment_score = get_final_score(genes, "L", 0.23, score, dataset) 
                                                    
                    if compartment_score < 0:
                        final_scores_with_compartments[name] = compartment_score
                    else:
                        final_scores_with_compartments[name] = score
                    #print(name, genes, comp, compartment_score, score)
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
                    if "_C" in rxnId:
                        final_scores_with_compartments[rxnId] = final_scores[rxnId]
    # Adding artificial scores to these reactions (CO2, water and O2 exchanges and complex IV"
    final_scores_with_compartments["UF03227_TCE"] = 6
    final_scores_with_compartments["UF00806_M"] = 6
    final_scores_with_compartments["UF03314_TCE"] = 6
    final_scores_with_compartments["UF03382_TCE"] = 6

    reaction_scores_path = os.path.join(results_path, fungi_id + "_reactionScores.npy")
    np.save(reaction_scores_path, final_scores_with_compartments)
    #print(f"DEBUG: Saved reactionScores file to {reaction_scores_path}")

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
            reaction = model.reactions.get_by_id(key)
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
            reaction = model.reactions.get_by_id(key)
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
    annotation.to_csv(os.path.join(results_path, fungi_id + "_annotations.csv"), index=False)

    print("Scores calculated")

    return final_scores_with_compartments
