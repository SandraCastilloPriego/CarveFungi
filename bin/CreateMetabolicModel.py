import argparse
import os
import re
import statistics
import sys


def load_runtime_dependencies():
    """Load the scientific/runtime libraries needed by the metabolic model pipeline.

    This function keeps optional project dependencies lazy so the CLI can show help
    or fail well before importing heavy packages like COBRApy or libSBML.

    Raises:
        SystemExit: If a required dependency is missing from the Python environment.
    """
    global CarveMeFuncPool, pd, cobra, np, EggNogScoring, compartment_prediction, readSBMLFromFile, writeSBMLToFile

    try:
        import CarveMeFuncPool
        import pandas as pd
        import cobra
        import numpy as np
        import EggNogScoring
        import compartment_prediction
        import libsbml

        readSBMLFromFile = libsbml.readSBMLFromFile
        writeSBMLToFile = libsbml.writeSBMLToFile
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing required dependency for this pipeline: "
            f"{exc.name}. Install the project requirements and try again."
        ) from exc


def parse_cli_args():
    """Parse and validate the command-line arguments required by the pipeline.

    Returns:
        argparse.Namespace: A namespace containing the validated fungal identifier,
            EggNOG annotation file path, and FASTA file path.

    Raises:
        SystemExit: If required arguments are missing, empty, or point to non-existent files.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a fungal metabolic model from EggNOG annotations and FASTA-based "
            "compartment predictions."
        ),
        epilog=(
            "Example:\n"
            "  python bin/CreateModelEggNogPool.py "
            "Saccharomyces_cerevisiae_S288C "
            "data/annotations/functional/annotations.emapper.annotations "
            "data/annotations/fasta/Saccharomyces_cerevisiae_S288C.faa"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "fungi_id",
        help="Unique identifier for the fungus/model output (for example, Saccharomyces_cerevisiae_S288C).",
    )
    parser.add_argument(
        "eggnog_file",
        help="Path to the EggNOG annotation file generated for the organism.",
    )
    parser.add_argument(
        "fasta_file",
        help="Path to the FASTA file used for compartment prediction and gene mapping.",
    )
    args = parser.parse_args()

    if not args.fungi_id.strip():
        parser.error("The fungi_id argument cannot be empty.")
    if not os.path.isfile(args.eggnog_file):
        parser.error(f"EggNOG annotation file not found: {args.eggnog_file}")
    if not os.path.isfile(args.fasta_file):
        parser.error(f"FASTA file not found: {args.fasta_file}")

    return args



def add_missing_reactions(name, model, universalModel):
    """Ensure the reconstructed model contains essential exchange and biomass reactions.

    Some required reactions may be absent when annotations are sparse. This helper adds
    the most important universal reactions and rebuilds a biomass objective tailored to
    the fungal model.

    Args:
        name: Fungus identifier used to decide whether a complex I reaction should be removed.
        model: The reconstructed COBRA model that may be missing key reactions.
        universalModel: The universal reference model used as the source of required reactions.

    Returns:
        cobra.Model: The updated model with missing reactions and a corrected biomass reaction.
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
    biomass = universalModel.reactions.get_by_id("BIOMASS_C")
    model.remove_reactions(["BIOMASS_C"])
    reaction = cobra.Reaction('BIOMASS_C')
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
        if "lipid_c" in met.id:
            new_biomass[met]=-0.075
        if "protein_1g_c" in met.id:
            new_biomass[met]=-0.36
        if "protein_1g_m" in met.id:
            new_biomass[met]=-0.04
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
    model.objective="BIOMASS_C"
    
    #Removing Complex I
    if "accharomyces" in name or "andida_glabrata" in name or "lyveromyces" in name or "shbya" in name or "chizosaccharomyces" in name:
        if model.reactions.has_id("UF02391_M"):
            model.remove_reactions(["UF02391_M"])
    else:
        if model.reactions.has_id("UF02391_M")==False:
            model.add_reactions([universalModel.reactions.get_by_id("UF02391_M")])

    return model


def combine_models(fungi_id,model_constrained, model_open_bounds, universal_model):
    """Merge the constrained and open-bound model ensembles into a single model.

    This combines the reactions that are unique to the open-bound reconstruction with the
    constrained reconstruction and marks each reaction as coming from the constrained or
    open model using the reaction name suffix.

    Args:
        fungi_id: Identifier for the fungal model being reconstructed.
        model_constrained: The model generated with constrained exchange bounds.
        model_open_bounds: The model generated with open exchange bounds.
        universal_model: The universal reference model used to fill missing reactions.

    Returns:
        cobra.Model: A merged model that includes both model variants.
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
    """Attach gene-protein-reaction rules to reactions based on EC-to-gene mappings.

    The function matches each EC entry to the universal reaction annotation table and
    appends matching gene IDs to each reaction rule using an OR relationship.

    Args:
        model: The COBRA model whose reactions will receive gene rules.
        ecs_to_genes_dictionary: EC number to gene/compartment mappings generated during scoring.
        universal: DataFrame containing the universal reaction annotation table.

    Returns:
        cobra.Model: The model with updated gene reaction rules.
    """
    ecs_to_genes_fixed = {}
    for key, item in ecs_to_genes_dictionary.items():
        genes = []
        for gen in item:
            if '-' in str(gen[0]):
                genes.append([gen[0].replace('-', ''),gen[1]])
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
                    for gene_comp in item:
                        compartments = gene_comp[1]                     
                        print(id_reaction, gene_comp)
                        for rxn in model.reactions:
                            iscompartment=False
                            for comp in compartments: 
                                if comp == 'O':
                                    comp = 'C'                         
                                if comp in rxn.id.split("_")[1]:
                                    iscompartment=True
                            if name in rxn.id and iscompartment and gene_comp[0] not in rxn.gene_reaction_rule:
                                rxn.gene_reaction_rule = rxn.gene_reaction_rule + " or " + gene_comp[0] if rxn.gene_reaction_rule != "" else gene_comp[0]
    return model


def annotate(fungi_id, doc, universal, compounds_mapping, objective, results_path="../results"):
    """Annotate the SBML model with compartment information, KEGG/ChEBI/PubChem metadata, and SBO terms.

    Args:
        fungi_id: Model identifier used in the metadata and output filename.
        doc: libSBML document representing the model to annotate.
        universal: DataFrame of universal reaction metadata used for reaction annotations.
        compounds_mapping: DataFrame mapping metabolite IDs to external identifiers.
        objective: Objective label to include in the output file name.
        results_path: Directory where the annotated SBML output is written.

    Returns:
        None. An annotated SBML file is written to disk.
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
        sel = compounds_mapping[compounds_mapping['IDs'] == str(spID)]
        if len(sel) == 0:
            spID = spID[0:len(spID) - 2]
            spID = spID + '_x'
            sel = compounds_mapping[compounds_mapping['IDs'] == str(spID)]
            if len(sel) == 0:
                spID = spID[0:len(spID) - 2]
                spID = spID + '_r'
                sel = compounds_mapping[compounds_mapping['IDs'] == str(spID)]
                if len(sel) == 0:
                    spID = spID[0:len(spID) - 2]
                    spID = spID + '_m'
                    sel = compounds_mapping[compounds_mapping['IDs'] == str(spID)]
                    if len(sel) == 0:
                        spID = spID[0:len(spID) - 2]
                        spID = spID + '_e'
                        sel = compounds_mapping[compounds_mapping['IDs'] == str(spID)]
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

    writeSBMLToFile(doc, os.path.join(results_path, fungi_id + '-' + objective + '.sbml'))
    os.system('sed -i \'s/<fbc:geneProduct /<fbc:geneProduct sboTerm="SBO:0000243" /g\' %s/%s-%s.sbml' % (results_path,
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


def carveFungi(open_bounds, fungi_id, universal_model_path, reaction_scores, ecs_genes_dict,universal, compounds_map, results_path):
    """Run the CarveMe reconstruction and save the per-model SBML outputs.

    Args:
        open_bounds: If True, use open exchange bounds; if False, use constrained bounds.
        fungi_id: Identifier for the fungal organism.
        universal_model_path: Path to the universal metabolic model used for carving.
        reaction_scores: Reaction-level scoring table used to prune the universal model.
        ecs_genes_dict: EC-to-gene mappings used to add gene rules to reconstructed reactions.
        universal: Universal reaction metadata table.
        compounds_map: Mapping of metabolite IDs to annotation identifiers.
        results_path: Directory where intermediate model files are saved.

    Returns:
        tuple: A three-item tuple containing the number of reconstructed models, the mean
            objective value, and the reconstructed model list.
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
        universal_model.reactions.get_by_id('UF02997_E').lower_bound = -100
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
            cobra.io.write_sbml_model(modelCreated, os.path.join(results_path, fungi_id + str(reconstructed_counter) + 'M.sbml'))

            doc = readSBMLFromFile(os.path.join(results_path, fungi_id + str(reconstructed_counter) + 'M.sbml'))
            annotate(fungi_id + "-M", doc, universal, compounds_map, str(objective[reconstructed_counter]), results_path)
            obj = obj + objective[reconstructed_counter]
            reconstructed_counter = reconstructed_counter + 1
        obj = obj / len(objective)

    else:
        reconstructed_counter = 0
        obj = 0
        for modelCreated in reconstructed_models:
            model_gen_annotations = add_gene_annotation(modelCreated, ecs_genes_dict, universal)
            modelCreated = add_missing_reactions(fungi_id, model_gen_annotations, universal_model)
            cobra.io.write_sbml_model(modelCreated, os.path.join(results_path, fungi_id + str(reconstructed_counter) + 'O.sbml'))

            doc = readSBMLFromFile(os.path.join(results_path, fungi_id + str(reconstructed_counter) + 'O.sbml'),)
            annotate(fungi_id + "-O",doc, universal, compounds_map, str(objective[reconstructed_counter]), results_path)
            obj = obj + objective[reconstructed_counter]
            reconstructed_counter = reconstructed_counter + 1
        obj = obj / len(objective)

    return reconstructed_counter, obj, reconstructed_models


def carveFungi_pipeline(fungi_id="test", eggnog_file=None, fasta_file = None):
    """Execute the full fungal model reconstruction workflow.

    The pipeline performs the following stages:
        1. predict protein compartment localization from the FASTA file,
        2. score reactions based on EggNOG annotations,
        3. run CarveMe using constrained and open-bound models,
        4. merge the reconstructed ensembles,
        5. write the final annotated model to the results directory.

    Args:
        fungi_id: Unique identifier for the fungal genome or model.
        eggnog_file: Path to the EggNOG annotation file.
        fasta_file: Path to the input FASTA file used for compartment prediction.

    Returns:
        None. The final SBML model and intermediate artifacts are written to the results folder.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    load_runtime_dependencies()

    # Construct the path to the data directory
    data_path = os.path.join(script_dir, '..', 'data')
    data_path = os.path.normpath(data_path)

    results_path = os.path.join(script_dir, '..', 'results')
    results_path = os.path.normpath(results_path)

    universal_model_path = os.path.join(data_path,'reactionDatabase/bigModelv3.8b.sbml')
    universal_csv_path = os.path.join(data_path,'reactionDatabase/universal_v3.8.csv')
    compounds_information_path = os.path.join(data_path,'reactionDatabase/NameMappingFinalv5.csv')
    
    print("Doing Compartment Prediction")

    if fasta_file !=None:
        compartment_file=compartment_prediction.run_compartment_predictions(os.path.join(data_path,'annotations/localization/'), fasta_file)
        if compartment_file is None:
            print("Error running compartment prediction. Exiting pipeline.")
            return
    elif os.path.isfile(data_path+'annotations/localization/' + fungi_id + ".loc_pred"):
        compartment_file = data_path+'annotations/localization/' + fungi_id + ".loc_pred"
    else:
        print("No FASTA or compartment prediction file provided, skipping reconstruction.")
        return
    
    print("Compartment Prediction Done")

    # Getting the reaction scores
    ec_to_score3, ecs_genes_dict, = EggNogScoring.get_gene_to_score_dict(
        fungi_id, eggnog_file, compartment_file)
    
    universal_data = EggNogScoring.assign_scores(ec_to_score3)
    final_scores = EggNogScoring.scoring_compartments(compartment_file, universal_data, ecs_genes_dict, fungi_id)
    

    universal = pd.read_csv(universal_csv_path)
    compounds_map = pd.read_csv(compounds_information_path)

    # Calling CarveMe
    model_files = os.listdir(results_path)
    M_models=[]
    for m in model_files:
        if fungi_id + "-M-" in m:
            M_models.append(m)

    if len(M_models)==0:
        number_of_models, objective, reconstructed_models = carveFungi(False, fungi_id, universal_model_path, final_scores, ecs_genes_dict,universal, compounds_map, results_path)
        if number_of_models == 0:
            print("Error: No valid models were reconstructed by CarveMe; check reaction scores, EC coverage, or model constraints.")
            return
    else:
        number_of_models=len(M_models)
        objective=[]
        reconstructed_models = []
        for m in M_models:
            objective.append(re.search("-M-\d+.\d+", m).group(0).replace("-M-",""))
            reconstructed_models.append(cobra.io.read_sbml_model(os.path.join(results_path,m)))
        obj = np.array(objective, dtype=float)
        objective= np.mean(obj)


    number_of_open_models, objective_open, reconstructed_models_open = carveFungi(True, fungi_id, universal_model_path, final_scores, ecs_genes_dict,universal, compounds_map,results_path)
    if number_of_open_models == 0:
        print("Error: No valid models were reconstructed by CarveMe; check reaction scores, EC coverage, or model constraints..")
        return

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

    try:
        assert final_model is not None, "Error: Merge step skipped because neither the constrained nor open model set was valid."
    except AssertionError as e:
        print(e)
        return
    
    for reaction in final_model.reactions:
        reaction.lower_bound = universal_model.reactions.get_by_id(reaction.id).lower_bound
        reaction.upper_bound = universal_model.reactions.get_by_id(reaction.id).upper_bound

    fix_ensemble_model(final_model)
    cobra.io.write_sbml_model(final_model, os.path.join(results_path, fungi_id + '-genes.sbml'))
    doc = readSBMLFromFile(os.path.join(results_path, fungi_id + '-genes.sbml'))
    annotate(fungi_id , doc, universal, compounds_map, "-" + str(objective) + "-" + str(objective_open), results_path)

    os.remove(os.path.join(results_path, fungi_id + '-genes.sbml'))
    for m in range(number_of_models):
        os.remove(os.path.join(results_path, fungi_id + str(m) + 'M.sbml'))
    for m in range(number_of_open_models):
        os.remove(os.path.join(results_path, fungi_id + str(m) + 'O.sbml'))

    print("Done", objective, objective_open)


def __main__():
    """Command-line entry point for running the full reconstruction workflow.

    The script expects three positional arguments:
        fungi_id, eggnog_file, and fasta_file.
    """
    args = parse_cli_args()
    load_runtime_dependencies()

    print("Starting the pipeline")
    carveFungi_pipeline(args.fungi_id, args.eggnog_file, args.fasta_file)


if __name__ == "__main__":
    __main__()