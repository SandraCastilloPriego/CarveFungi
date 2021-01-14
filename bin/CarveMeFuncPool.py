'''
  CarveMe algorithm. Code extracted from: https://github.com/cdanielmachado/carveme								
  D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018..
  DOI: https://doi.org/10.1093/nar/gky537							
  Developed at the European Molecular Biology Laboratory (2017-2018). Released under an Apache License. 
   
'''

from framed import Parameter
from framed.cobra.ensemble import EnsembleModel, save_ensemble
from framed.io.sbml import parse_gpr_rule, save_cbmodel
from framed.model.transformation import disconnected_metabolites
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from framed import FBA
from collections import Counter
import operator


def minmax_reduction(model, scores, min_growth=0.1, min_atpm=0.1, eps=1e-5, bigM=1e3, default_score=-1.0,
                     uptake_score=1, soft_score=1.0, soft_constraints=[], hard_constraints=None,
                     ref_reactions=[], ref_score=0.0, solver=None, debug_output=None, feast=1e-6, opti=1e-5):
    """ Apply minmax reduction algorithm (MILP).
    Computes a binary reaction vector that optimizes the agreement with reaction scores (maximizes positive scores,
    and minimizes negative scores). It generates a fully connected reaction network (i.e. all reactions must be able
    to carry some flux).
    Args:
        model (CBModel): universal model
        scores (dict): reaction scores
        min_growth (float): minimal growth constraint
        min_atpm (float): minimal maintenance ATP constraint
        eps (float): minimal flux required to consider leaving the reaction in the model
        bigM (float): maximal reaction flux
        default_score (float): penalty score for reactions without an annotation score (default: -1.0).
        uptake_score (float): penalty score for using uptake reactions (default: 0.0).
        soft_score (float): score for soft constraints (default: 1.0)
        soft_constraints (dict): dictionary from reaction id to expected flux direction (-1, 1, 0)
        hard_constraints (dict): dictionary of flux bounds
        solver (Solver): solver instance (optional)
    Returns:
        Solution: optimization result
    """

    if not solver:
        solver = solver_instance(model)

    objective = {}
    scores = scores.copy()
    reactions = []
    for r_id in model.reactions:
        if r_id not in reactions and not r_id.endswith('_E'):
            reactions.append(r_id)

    if soft_constraints:
        reactions += [r_id for r_id in soft_constraints if r_id not in reactions]
    else:
        soft_constraints = {}

    if hard_constraints:
        solver.set_bounds(hard_constraints)

    # R_UF01847_CE is the ATP maintenance reaction
    # if the default score is lower than 0 set all the reactions to the default score except for the exchange reactions
    # and the ATP maintenance reaction

    if default_score != 0:
        for r_id in model.reactions:
            if r_id not in reactions and r_id not in ref_reactions and not r_id.endswith(
                    '_E') and r_id != 'R_UF01847_CE':
                scores[r_id] = default_score
                reactions.append(r_id)

    if ref_score != 0:
        for r_id in ref_reactions:
            if r_id not in reactions and r_id != 'R_UF01847_CE':
                scores[r_id] = ref_score
                reactions.append(r_id)

    if not hasattr(solver, '_carveme_flag'):
        solver._carveme_flag = True

        solver.add_constraint('min_growth', {'R_BIOMASS': 1}, '>', min_growth, update_problem=False)
        solver.add_constraint('min_atpm', {'R_UF01847_CE': 1}, '>', min_atpm, update_problem=False)

        solver.neg_vars = []
        solver.pos_vars = []

        for r_id in reactions:
            if model.reactions[r_id].lb is None or model.reactions[r_id].lb < 0:
                y_r = 'yr_' + r_id
                solver.add_variable(y_r, 0, 1, vartype=VarType.BINARY, update_problem=False)
                solver.neg_vars.append(y_r)
            if model.reactions[r_id].ub is None or model.reactions[r_id].ub > 0:
                y_f = 'yf_' + r_id
                solver.add_variable(y_f, 0, 1, vartype=VarType.BINARY, update_problem=False)
                solver.pos_vars.append(y_f)

        if uptake_score != 0:
            for r_id in model.reactions:
                if r_id.endswith('_E'):
                    solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update_problem=False)

        solver.update()

        for r_id in reactions:
            y_r, y_f = 'yr_' + r_id, 'yf_' + r_id
            if y_r in solver.neg_vars and y_f in solver.pos_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_f: -eps, y_r: bigM}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_f: -bigM, y_r: eps}, '<', 0, update_problem=False)
                solver.add_constraint('rev_' + r_id, {y_f: 1, y_r: 1}, '<', 1, update_problem=False)
            elif y_f in solver.pos_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_f: -eps}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_f: -bigM}, '<', 0, update_problem=False)
            elif y_r in solver.neg_vars:
                solver.add_constraint('lb_' + r_id, {r_id: 1, y_r: bigM}, '>', 0, update_problem=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, y_r: eps}, '<', 0, update_problem=False)

        if uptake_score != 0:
            for r_id in model.reactions:
                if r_id.endswith('_E'):
                    solver.add_constraint('lb_' + r_id, {r_id: 1, 'y_' + r_id: bigM}, '>', 0, update_problem=False)

        solver.update()

    for r_id in reactions:
        y_r, y_f = 'yr_' + r_id, 'yf_' + r_id

        if r_id in soft_constraints:
            sign = soft_constraints[r_id]

            if sign > 0:
                w_f, w_r = soft_score, 0
            elif sign < 0:
                w_f, w_r = 0, soft_score
            else:
                w_f, w_r = -soft_score, -soft_score

        if y_f in solver.pos_vars:
            if r_id in soft_constraints:
                objective[y_f] = w_f
            elif ref_score != 0 and r_id in ref_reactions:
                objective[y_f] = 2 * scores[r_id] + ref_score
            else:
                objective[y_f] = scores[r_id]

        if y_r in solver.neg_vars:
            if r_id in soft_constraints:
                objective[y_r] = w_r
            elif ref_score != 0 and r_id in ref_reactions:
                objective[y_r] = 2 * scores[r_id] + ref_score
            else:
                objective[y_r] = scores[r_id]

    if uptake_score != 0:
        for r_id in model.reactions:
            if r_id.endswith('_E') and r_id not in soft_constraints:
                objective['y_' + r_id] = uptake_score

    if debug_output:
        solver.write_to_file(debug_output + "_milp_problem.lp")

    solver.set_parameter(Parameter.INT_FEASIBILITY_TOL, feast)
    solver.set_parameter(Parameter.OPTIMALITY_TOL, opti)

    solutions = solver.solve(linear=objective, minimize=False, get_values=True, pool_size=50, pool_gap=0)

    return solutions


def carve_model(model, reaction_scores, outputfile=None, min_growth=0.1, min_atpm=0.1, flavor=None, inplace=False,
                default_score=-0.1, uptake_score=1.0, soft_score=1.0,
                soft_constraints=None, hard_constraints=None,
                init_env=None, debug_output=None, eps=1e-5, feast=1e-6, opti=1e-6):
    """ Reconstruct a metabolic model using the CarveMe approach.
    Args:
        model (CBModel): universal model
        reaction_scores (pandas.DataFrame): reaction scores
        outputfile (str): write model to SBML file (optional)
        flavor (str): SBML flavor ('cobra' or 'fbc2', optional)
        inplace (bool): Change model in place (default: True)
        default_score (float): penalty for non-annotated intracellular reactions (default: -1.0)
        uptake_score (float): penalty for utilization of extracellular compounds (default: 0.0)
        soft_score (float): score for soft constraints (default: 1.0)
        soft_constraints (dict): dictionary from reaction id to expected flux direction (-1, 1, 0)
        hard_constraints (dict): dictionary of flux bounds
        init_env (Environment): initialize final model with given Environment (optional)
    Returns:
        CBModel: reconstructed model
    """

    scores = reaction_scores

    if default_score is None:
        default_score = -1.0

    if uptake_score is None:
        uptake_score = 0.0

    if soft_score is None:
        soft_score = 1.0

    if soft_constraints:
        not_in_model = set(soft_constraints) - set(model.reactions)
        if not_in_model:
            soft_constraints = {r_id: val for r_id, val in soft_constraints.items() if r_id in model.reactions}
            warnings.warn("Soft constraints contain reactions not in the model:\n" + "\n".join(not_in_model))

    if hard_constraints:
        not_in_model = set(hard_constraints) - set(model.reactions)
        if not_in_model:
            hard_constraints = {r_id: (lb, ub) for r_id, (lb, ub) in hard_constraints.items() if
                                r_id in model.reactions}
            warnings.warn("Hard constraints contain reactions not in the model:\n" + "\n".join(not_in_model))
    reconstructed_model = model.copy()
    solutions = minmax_reduction(reconstructed_model, scores, eps=eps, default_score=default_score, uptake_score=uptake_score,
                                 soft_score=soft_score, soft_constraints=soft_constraints, min_growth=min_growth,
                                 min_atpm=min_atpm,
                                 hard_constraints=hard_constraints, debug_output=debug_output, feast=feast, opti=opti)

    solutions = select_best_models(solutions)

    models = []
    i = 0
    objective = 0
    try:
        for sol in solutions:
            reconstructed_model = model.copy()
            print("Solution " + str(i))
            print(sol)
            i = i + 1
            inactive = inactive_reactions(reconstructed_model, sol, reaction_scores)
            reconstructed_model.remove_reactions(inactive)
            del_metabolites = disconnected_metabolites(reconstructed_model)
            reconstructed_model.remove_metabolites(del_metabolites)
            print(FBA(reconstructed_model, objective={'R_BIOMASS': 1}))
            models.append(reconstructed_model)
            objective = objective + sol.fobj
        objective = objective / len(solutions)
    except:
        reconstructed_model = model.copy()
        print(solutions)
        inactive = inactive_reactions(reconstructed_model, solutions, reaction_scores)
        reconstructed_model.remove_reactions(inactive)
        del_metabolites = disconnected_metabolites(reconstructed_model)
        reconstructed_model.remove_metabolites(del_metabolites)
        print(FBA(reconstructed_model, objective={'R_BIOMASS': 1}))
        models.append(reconstructed_model)
        objective = solutions.fobj

    return objective, models


def select_best_models(solutions):
    best_candidates = {}
    for sol in solutions:
        best_candidates[sol] = sol.fobj
    best_solutions = dict(sorted(best_candidates.items(), key=operator.itemgetter(1), reverse=True)[0:5])
    return best_solutions


def inactive_reactions(model, solution, score):
    inactive = []

    internal = [r_id for r_id in model.reactions if not r_id.endswith('_E')]
    external = [r_id for r_id in model.reactions if r_id.endswith('_E')]

    variable = []
    sol = []
    for r_id in internal:
        threshold = 0.5
        variable.append(r_id)
        sol.append(solution.values[r_id])
        if solution.values.get('yf_' + r_id, 0) < threshold and solution.values.get('yr_' + r_id, 0) < threshold:
            inactive.append(r_id)
        variable.append('yf_' + r_id)
        sol.append(solution.values.get('yf_' + r_id, 0))
        variable.append('yr_' + r_id)
        sol.append(solution.values.get('yr_' + r_id, 0))
    m_r_lookup = model.metabolite_reaction_lookup()
    inactive_ext = []

    for r_id in external:
        try:
            m_id = model.reactions[r_id].get_substrates()[0]
            neighbors = m_r_lookup[m_id]
            if len(set(neighbors) - set(inactive)) == 1:
                inactive_ext.append(r_id)
        except:
            print(r_id, m_id, "is there something wrong?")
    return inactive + inactive_ext
