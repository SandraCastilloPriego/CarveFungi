import CarveMeFuncPool
import pandas as pd
import cobra
import re
from libsbml import *
import EggNogScoring
import statistics
import sys
import os
from cobra import Model, Reaction, Metabolite


def add_missing_reactions(name, model, universalModel):
    """
    Adds some reactions that should be in all the models in case they are not added during the reconstruction
    due to the lack of annotations.
    :param model: Cobra model.
    :param universalModel: Universal cobra model.
    :return: the Cobra model with the reactions added.
    """
    # Adding CO2 exchange
    if model.reactions.has_id("UF03227_E") == False or model.reactions.has_id("UF03227_TCE") == False:
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF03227_TCE")])
        except:
            print("Not able to add CO2 exchange")
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF03227_E")])
        except:
            print("Not able to add CO2 exchange")
    # Adding Ammonium exchange
    if model.reactions.has_id("UF03376_E") == False or model.reactions.has_id("UF03376_TCE") == False:
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF03376_TCE")])
        except:
            print("Not able to add Ammonium exchange")
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF03376_E")])
        except:
            print("Not able to add Ammonium exchange")
    # Adding complex IV
    if not model.reactions.has_id("UF00806_M"):
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF00806_M")])
        except:
            print("Not able to add complex IV")
    # Adding water exchange
    if not model.reactions.has_id("UF03314_TCE"):
        try:
            model.add_reaction(
                universalModel.reactions.get_by_id("UF03314_TCE"))
        except:
            print("Not able to add water exchange")
        try:
            model.add_reactions([universalModel.reactions.get_by_id("UF03314_E")])
        except:
            print("Not able to add water exchange")


    # Fixing Biomass
    biomass = universalModel.reactions.get_by_id("BIOMASS")
    model.remove_reactions(["BIOMASS"])
    reaction = Reaction('BIOMASS')
    reaction.name = 'BIOMASS - artificial'
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.
    new_biomass = {}

    for met in biomass.metabolites:
        if "13BDglcn_1g_c" in met.id:
            new_biomass[met]=-0.253
        if "ash_1g_c" in met.id:
            new_biomass[met]=-0.048
        if "atp_c" in met.id:
            new_biomass[met]=-60.0
        if "chitin_1g_c" in met.id:
            new_biomass[met]=-0.01
        if "dna_1g_c" in met.id:
            new_biomass[met]=-0.005
        if "glycogen_1g_c" in met.id:
            new_biomass[met]=-0.1
        if "h2o_c" in met.id:
            new_biomass[met]=-60
        if "mannan_1g_c" in met.id:
            new_biomass[met]=-0.07
        if "phospholipid_1g_c" in met.id:
            new_biomass[met]=-0.04
        if "protein_1g_c" in met.id:
            new_biomass[met]=-0.4
        if "rna_1g_c" in met.id:
            new_biomass[met]=-0.058
        if "sterol_1g_r" in met.id:
            new_biomass[met]=-0.001
        if "trehalose_1g_c" in met.id:
            new_biomass[met]=-0.015
        if "ficytc_c" in met.id:
            new_biomass[met]=-0.0001
        if "focytc_c" in met.id:
            new_biomass[met]=-0.0001
        if "retn_c" in met.id:
            new_biomass[met]=-0.0001
        if "thmtp_c" in met.id:
            new_biomass[met]=-0.0001
        if "ribflv_c" in met.id:
            new_biomass[met]=-0.0001
        if "btn_c" in met.id:
            new_biomass[met]=-0.0001
        if "q6_m" in met.id:
            new_biomass[met]=-0.0001
        if "coa_c" in met.id:
            new_biomass[met]=-0.0001
        if "thf_c" in met.id:
            new_biomass[met]=-0.0001
        if "adp_c" in met.id:
            new_biomass[met]=60.0
        if "biomass_1g_c" in met.id:
            new_biomass[met]=1
        if "h_c" in met.id:
            new_biomass[met]=60.0
        if "pi_c" in met.id:
            new_biomass[met]=60.0
    reaction.add_metabolites(new_biomass)
    model.add_reactions([reaction])
    model.objective="BIOMASS"
    '''
    #Adding Complex I
    if "accharomyces" in name or "andida_glabrata" in name or "lyveromyces" in name or "shbya" in name or "chizosaccharomyces" in name:
        if model.reactions.has_id("UF02391_M"):
            model.remove_reactions(["UF02391_M"])
    else:
        if model.reactions.has_id("UF02391_M")==False:
            model.add_reactions([universalModel.reactions.get_by_id("UF02391_M")])'''

    return model


def combine_models(fungi_id,model_constrained, model_open_bounds, universal_model):
    """
    Combines the model assemblies coming from the two CarveMe simulations.
    :param model_constrained: Cobra model obtained from the CarveMe simulation using the constrained universal model.
    :param model_open_bounds: Cobra model obtained from the CarveMe simulation using the universal model without
        constrains.
    :param universal_model:
    :return: The combined cobra model.
    """
    reaction_differences = []
    model_open_bounds = add_missing_reactions(fungi_id, model_open_bounds, universal_model)

    for reaction in model_open_bounds.reactions:
        if reaction not in model_constrained.reactions:
            reaction_differences.append(reaction)

    reaction_differences = list(set(reaction_differences))
    model_constrained.add_reactions(reaction_differences)

    for reaction in model_constrained.reactions:
        if model_open_bounds.reactions.has_id(reaction.id):
            reaction.name = reaction.name + "_" + str(1)
        else:
            reaction.name = reaction.name + "_" + str(0)

    return model_constrained


def add_gene_annotation(model, ecs_to_genes_dictionary, universal):
    """
    Adds the gene annotation based on the EC to genes dictionary created during the reaction scoring. The rules
    are unknown so in all the cases the genes are added with the rule "OR". This needs to be changed when more
    information about the models is collected.
    :param model: Cobra model.
    :param ecs_to_genes_dictionary: dictionary containing the correspondence between ECs and genes.
    :param universal: Pandas dataFrame containing all the annotation information.
    :return: Cobra models with the genes-reaction association.
    """
    ecs_to_genes_fixed = {}
    for key, item in ecs_to_genes_dictionary.items():
        genes = []
        for gen in item:
            if '-' in str(gen):
                genes.append(gen.replace('-', ''))
            else:
                genes.append(gen)
        ecs_to_genes_fixed[key] = genes

    for key, item in ecs_to_genes_fixed.items():
        if len(item) > 0:
            sel = universal[universal["ECs"].str.contains(key)]
            for k, r in sel.iterrows():
                ECS=r["ECs"].split("|")
                goahead=False
                for EC in ECS:
                    if key == EC:
                        goahead=True
                if goahead:
                    id_reaction = r["IDs"]
                    name = id_reaction.split("_")[0]
                    for rxn in model.reactions:
                        if name in rxn.id and len(rxn.genes) < len(item):
                            rxn.gene_reaction_rule = (" or ".join(item))
    return model


def annotate(fungi_id, doc, universal, compounds_mapping, objective):
    """
    Adds the annotations to the sbml file.
    :param fungi_id: Fungi identifier.
    :param doc: Object of the sbml file document.
    :param universal: Pandas dataframe containing the reaction information.
    :param compounds_mapping: Pandas dataframe containing the compound information.
    :param objective: String with the objective information. It will be added to the file name.
    """
    model = doc.getModel()
    model.setId(fungi_id)
    model.setMetaId('meta_' + fungi_id)

    NOTES = """
    <body xmlns="http://www.w3.org/1999/xhtml">
    <p>
    Metabolic model of %s.
    This model was reconstructed with CarveMe for Fungi.
    </p>
    </body>
    """ % (fungi_id)

    model.setNotes(NOTES)

    for compound in model.getListOfSpecies():
        # print(compound.getId())
        # print(compound.getName())
        spID = compound.getId()
        spID = spID[2:len(spID) - 2]
        spID = spID + '_c'
        sel = compounds_mapping[compounds_mapping['id'] == str(spID)]
        if len(sel) == 0:
            spID = spID[0:len(spID) - 2]
            spID = spID + '_x'
            sel = compounds_mapping[compounds_mapping['id'] == str(spID)]
            if len(sel) == 0:
                spID = spID[0:len(spID) - 2]
                spID = spID + '_r'
                sel = compounds_mapping[compounds_mapping['id'] == str(spID)]
                if len(sel) == 0:
                    spID = spID[0:len(spID) - 2]
                    spID = spID + '_m'
                    sel = compounds_mapping[compounds_mapping['id'] == str(spID)]
                    if len(sel) == 0:
                        spID = spID[0:len(spID) - 2]
                        spID = spID + '_e'
                        sel = compounds_mapping[compounds_mapping['id'] == str(spID)]
                        if len(sel) == 0:
                            print(spID)
        compound.setMetaId(compound.getId())
        compound.unsetAnnotation()
        compound.setSBOTerm("SBO:0000247")
        if len(sel) > 0:
            kegg = sel.iloc[0]['Kegg']
            chebi = sel.iloc[0]['ChEBI']
            pubchem = str(sel.iloc[0]['PubChem'])
            pubchem = pubchem[0:pubchem.find('.')]
            if type(kegg) != float or type(chebi) != float or (type(pubchem) != float and len(pubchem) > 1):
                annotation = """<annotation xmlns:sbml="http://www.sbml.org/sbml/level3/version1/core">
            <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
                <rdf:Description rdf:about="#%s">
                     <bqbiol:is>
                         <rdf:Bag>""" % (spID)
                if type(kegg) != float:
                    annotation = annotation + """
                             <rdf:li rdf:resource="http://identifiers.org/kegg.compound/%s"/>  """ % (kegg)
                if type(chebi) != float:
                    annotation = annotation + """
                             <rdf:li rdf:resource="http://identifiers.org/chebi/%s"/>""" % (chebi)
                if type(pubchem) != float and len(pubchem) > 1:
                    annotation = annotation + """
                             <rdf:li rdf:resource="http://identifiers.org/pubchem.compound/%s"/>""" % (pubchem)

                annotation = annotation + """
                         </rdf:Bag>
                     </bqbiol:is>
                     <bqbiol:hasProperty>
                        <rdf:Bag>
                          <rdf:li rdf:resource="http://identifiers.org/sbo/SBO:0000247"/>
                        </rdf:Bag>
                      </bqbiol:hasProperty>
                </rdf:Description>
            </rdf:RDF>
        </annotation>"""
                compound.setAnnotation(annotation)

                # print(compound.getAnnotationString())
        else:
            print(spID)

    for reaction in model.getListOfReactions():
        # print(reaction.getId())
        # print(reaction.getName())
        spID = reaction.getId()
        if '_x' in spID[len(spID) - 2:len(spID)] or '_m' in spID[len(spID) - 2:len(spID)] or '_c' in spID[len(spID) - 2:len(spID)] or '_r' in spID[len(spID) - 2:len(spID)]:
            spID = spID[0:len(spID) - 2]
        sel = universal[universal['IDs'] == spID[2:len(spID)]]
        reaction.setMetaId(reaction.getId())
        reaction.unsetAnnotation()
        if len(sel) > 0:
            if 'transport' in sel.iloc[0]['ECs']:
                reaction.setSBOTerm("SBO:0000185")
            elif "exchange" in sel.iloc[0]['ECs']:
                reaction.setSBOTerm("SBO:0000627")
            elif 'artificial' in sel.iloc[0]['ECs']:
                reaction.setSBOTerm("SBO:0000629")
            else:
                reaction.setSBOTerm("SBO:0000176")

            if "UF03247_C" in reaction.getId():
                reaction.setSBOTerm("SBO:0000628")
            kegg = str(sel.iloc[0]['DB_IDs'])
            ec = sel.iloc[0]['ECs']
            if kegg.startswith('R') or '.' in ec:
                annotation = """<annotation xmlns:sbml="http://www.sbml.org/sbml/level3/version1/core">
                  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
                    <rdf:Description rdf:about="#%s">""" % (reaction.getId())
                if kegg.startswith('R'):
                    annotation = annotation + """  <bqbiol:isDescribedBy>
                            <rdf:Bag>
                              <rdf:li rdf:resource="http://identifiers.org/kegg.reaction/%s"/>
                            </rdf:Bag>
                          </bqbiol:isDescribedBy>""""""
                          """ % (kegg)
                if '.' in ec:
                    annotation = annotation + """<bqbiol:is>
                            <rdf:Bag>  """
                    for ecnumber in ec.split('|'):
                        annotation = annotation + """
                                  <rdf:li rdf:resource="http://identifiers.org/ec-code:%s"/>
                        """ % (ecnumber)

                    annotation = annotation + """        </rdf:Bag>
                      </bqbiol:is>"""

                annotation = annotation + """
                    </rdf:Description>
                  </rdf:RDF>
                </annotation>"""
                reaction.setAnnotation(annotation)
            reaction.getSBOTerm()
        else:
            if '_t' in reaction.id:
                reaction.setSBOTerm("SBO:0000655")

    writeSBMLToFile(doc, '../results/' + fungi_id + '-' + objective + '.sbml')
    os.system('sed -i \'s/<fbc:geneProduct /<fbc:geneProduct sboTerm="SBO:0000243" /g\' ../results/%s-%s.sbml' % (
        fungi_id, objective))


def fix_ensemble_model(model):
    """
    this function adjusts the number of zeros and ones in the reaction names of the model.
    :param model: Model to be fixed.
    """
    count = []
    ensemble_reactions = {}
    for reaction in model.reactions:
        name = reaction.name.split('_')
        ensemble = []
        for n in name[1:]:
            try:
                n = int(n)
                ensemble.append(n)
            except:
                print(reaction.id, n)
        # print(ensemble)
        count.append(len(ensemble))
        ensemble_reactions[reaction.id] = ensemble
    mode = statistics.mode(count)
    for reaction in model.reactions:
        name = reaction.name.split('_')
        ensemble = []
        for n in name[1:]:
            try:
                n = int(n)
                ensemble.append(n)
            except:
                print(reaction.id, n)
        if len(ensemble) < mode:
            new_name = name[0] + "_"
            for i in range(mode - len(ensemble)):
                new_name = new_name + '0_'
            new_name = new_name + '_'.join(str(e) for e in ensemble)
            reaction.name = new_name


def carveFungi(open_bounds, fungi_id, universal_model_path, reaction_scores, ecs_genes_dict,universal, compounds_map):
    """
    Calls the CarveMe algorithm for the model reconstruction.
    :param open_bounds: True/False depending on the exchange constrains of the universal model.
    :param fungi_id: Fungi id.
    :param universal_model_path: Path to the universal sbml file.
    :param reaction_scores: Scores for the reactions.
    :return: The number of models reconstructed and the average objective value of the models in the ensemble.
    """
    # Open Framed model
    universal_model = cobra.io.read_sbml_model(universal_model_path)

    # Set the bounds of the model
    for reaction in universal_model.reactions:
        if reaction.lower_bound < 0:
            reaction.lower_bound = -100
        if reaction.upper_bound > 0:
            reaction.upper_bound = 100

    if not open_bounds:
        # Glucose, ammonium, water, O2,...
        universal_model.reactions.get_by_id('UF03376_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF02549_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03382_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03474_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF02765_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03268_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03456_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03314_E').lower_bound = -100
        universal_model.reactions.get_by_id('UF03288_E').lower_bound = -10
    else:
        for reaction in universal_model.reactions:
            if reaction.id.endswith('_E'):
                reaction.lower_bound = -100
                reaction.upper_bound = 100

    # Carving
    objective, reconstructed_models = CarveMeFuncPool.carve_model(
        universal_model, reaction_scores, eps=1e-3, min_growth=0.1, min_atpm=0.1, feast=1e-7, opti=1e-7)


    # Save the models into files (to be able to work with cobra and update the sbml files)
    if not open_bounds:
        reconstructed_counter = 0
        obj=0
        for modelCreated in reconstructed_models:
            model_gen_annotations = add_gene_annotation(modelCreated, ecs_genes_dict, universal)
            modelCreated = add_missing_reactions(fungi_id, model_gen_annotations, universal_model)
            cobra.io.write_sbml_model(modelCreated,'../results/' + fungi_id + str(reconstructed_counter) + 'M.sbml')

            doc = readSBMLFromFile('../results/' + fungi_id + str(reconstructed_counter) + 'M.sbml')
            annotate(fungi_id+ "-M", doc, universal, compounds_map, str(objective[reconstructed_counter]))
            obj = obj + objective[reconstructed_counter]
            reconstructed_counter = reconstructed_counter + 1
        obj = obj / len(objective)

    else:
        reconstructed_counter = 0
        obj = 0
        for modelCreated in reconstructed_models:
            model_gen_annotations = add_gene_annotation(modelCreated, ecs_genes_dict, universal)
            modelCreated = add_missing_reactions(fungi_id, model_gen_annotations, universal_model)
            cobra.io.write_sbml_model(modelCreated, '../results/' + fungi_id + str(reconstructed_counter) + 'O.sbml')

            doc = readSBMLFromFile('../results/' + fungi_id + str(reconstructed_counter) + 'O.sbml')
            annotate(fungi_id + "-O",doc, universal, compounds_map, str(objective[reconstructed_counter]))
            obj = obj + objective[reconstructed_counter]
            reconstructed_counter = reconstructed_counter + 1
        obj = obj / len(objective)

    return reconstructed_counter, obj, reconstructed_models


def carveFungi_pipeline(fungi_id="test", eggnog_file=None,  compartment_file=None):
    """
    Pipeline for the model reconstruction.

    :param fungi_id: Name of the model file.
    :param eggnog_file: Path to the functional annotation file (from eggNog).
    :param compartment_file: Path to the protein localization prediction.

    """
    universal_model_path = '../data/reactionDatabase/bigModelv2.21b.sbml'
    universal_csv_path = '../data/reactionDatabase/universal_v2.21.csv'
    compounds_information_path = '../data/reactionDatabase/NameMappingFinalv3.csv'

    if not os.path.isfile(compartment_file):
        sys.exit("Compartment file missing")

    # Getting the reaction scores
    ec_to_score3, ecs_genes_dict, = EggNogScoring.get_gene_to_score_dict(
        fungi_id, eggnog_file)
    if not os.path.isfile("../results/" + fungi_id + "_annotations.csv"):
        universal_data = EggNogScoring.assign_scores(ec_to_score3)
        final_scores = EggNogScoring.scoring_compartments(compartment_file, universal_data, ecs_genes_dict, fungi_id)
    else:
        scores = pd.read_csv("../results/" + fungi_id + "_annotations.csv")
        final_scores = {}
        for key, row in scores.iterrows():
            final_scores[row['IDs']] = row['Score']
    universal = pd.read_csv(universal_csv_path)
    compounds_map = pd.read_csv(compounds_information_path)

    # Calling CarveMe
    model_files = os.listdir("../results/")
    M_models=[]
    for m in model_files:
        if fungi_id + "-M-" in m:
            M_models.append(m)

    if len(M_models)==0:
        number_of_models, objective, reconstructed_models = carveFungi(False, fungi_id, universal_model_path, final_scores, ecs_genes_dict,universal, compounds_map)
    else:
        number_of_models=len(M_models)
        objective=[]
        reconstructed_models = []
        for m in M_models:
            objective.append(re.search("-M-\d+.\d+", m).group(0).replace("-M-",""))
            reconstructed_models.append(cobra.io.read_sbml_model("../results/"+m))

        objective= np.mean(np.array(objective))

    number_of_open_models, objective_open, reconstructed_models_open = carveFungi(True, fungi_id, universal_model_path, final_scores,ecs_genes_dict,universal, compounds_map)

    # Merging the models and adding the annotations to the sbml files.
    universal_model = cobra.io.read_sbml_model(universal_model_path)

    final_model = None

    for m1 in reconstructed_models:
        for m2 in reconstructed_models_open:
            model_merged = m1.merge(m2, inplace=False)
            if final_model is None:
                final_model = model_merged
                for reaction in final_model.reactions:
                    reaction.name = reaction.name + "_" + str(1)
            else:
                final_model = combine_models(fungi_id, final_model, model_merged, universal_model)

    for reaction in final_model.reactions:
        reaction.lower_bound = universal_model.reactions.get_by_id(reaction.id).lower_bound
        reaction.upper_bound = universal_model.reactions.get_by_id(reaction.id).upper_bound

    fix_ensemble_model(final_model)
    cobra.io.write_sbml_model(final_model, "../results/" + fungi_id + '-genes.sbml')
    doc = readSBMLFromFile('../results/' + fungi_id + '-genes.sbml')
    annotate(fungi_id , doc, universal, compounds_map, "-" + str(objective) + "-" + str(objective_open))

    os.remove('../results/'+fungi_id+ '-genes.sbml')
    for m in range(number_of_models):
        os.remove('../results/' + fungi_id + str(m) + 'M.sbml')
    for m in range(number_of_open_models):
        os.remove('../results/' + fungi_id + str(m) + 'O.sbml')

    print("Done", objective, objective_open)


def __main__():
    """
    The pipeline can be called from the command line. Example:
    data/annotations/functional/Saccharomyces_cerevisiae_S288C_GCF_000146045.2_NEWprotein.faa.emapper.annotations \
    data/annotations/fasta/Saccharomyces_cerevisiae_S288C_GCF_000146045.2_NEWprotein.faa \
    data/annotations/localization/Saccharomyces_cerevisiae_S288C_GCF_000146045.2_NEWprotein.faa.npy
    """
    fungi_id = sys.argv[1]
    eggnog_file = sys.argv[3]
    sequence_file = sys.argv[4]
    compartment_file = sys.argv[5]
    carveFungi_pipeline(fungi_id, eggnog_file, sequence_file, compartment_file)
