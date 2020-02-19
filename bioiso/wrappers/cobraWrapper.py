from cobra import io, Reaction
from cobra.flux_analysis import single_reaction_deletion, double_reaction_deletion
from bioiso.utils.bioisoUtils import NodeCache
import numpy as np
import random
import warnings
warnings.filterwarnings("ignore")


def load(file_name):
    try:
        with io.read_sbml_model(file_name) as model:
            return model

    except OSError:
        raise OSError("The directory or file name {} could not be found. Alternatively, this is not a valid SBML model. "
              "Consult COBRApy I/O for the supported formats".format(str(file_name)))

    except:
        print("Unexpected error")
        raise

def get_reaction(model, reaction_id):
    return model.reactions.get_by_id(reaction_id)

def has_reaction(model, reaction_id):
    return reaction_id in model.reactions

def get_metabolite(model, metabolite_id):
    return model.metabolites.get_by_id(metabolite_id)

def set_objective_function(model, reaction_id):
    reaction = get_reaction(model, reaction_id)
    model.objective = reaction

def set_solver(model, solver):

    try:
        model.solver = solver

    except:
        print('Using COBRApy default solver, namely glpk')

def objective(model):
    print(model.objective)

def get_reactants(model, reaction_id):
    return get_reaction(model, reaction_id).reactants

def get_products(model, reaction_id):
    return get_reaction(model, reaction_id).products

def list_reactants_names(model, reaction_id):
    reacs = get_reaction(model, reaction_id).reactants
    return [reac.name for reac in reacs]

def list_products_names(model, reaction_id):
    prods = get_reaction(model, reaction_id).products
    return [prod.name for prod in prods]

def get_reactions(model, metabolite_id):
    return get_metabolite(model, metabolite_id).reactions

def get_biomass_macromolecules(model,biomass_reaction_id):
    return get_reactants(model, biomass_reaction_id)

def create_unbalenced_reaction(model,metabolite_id, bounds=(-999999, 999999)):

    # bounds = (0, 999999)
    # demand - unbalanced network reaction that only allows the accumulation of a compound

    # bounds = (-999999, 0)
    # demand - unbalanced network reaction that only allows the uptake of a compound

    # bounds = (-999999, 999999)
    # sink - unbalanced network reaction that both provides or consumes the network with metabolites

    if bounds[0] >= 0.0:
        reaction_name = "Demand_" + metabolite_id
    else:
        reaction_name = "Sink_" + metabolite_id

    if has_reaction(model, reaction_name):
        return reaction_name

    else:

        model.add_reaction(Reaction(reaction_name))
        get_reaction(model,reaction_name).add_metabolites({get_metabolite(model, metabolite_id): -1})
        get_reaction(model,reaction_name).bounds = bounds

        return reaction_name

def get_reactions_by_role(model, metabolite_id, isReactant, previous_reactions_list):

    reactions = get_reactions(model, metabolite_id)

    reactions_list = []
    other_reactions_list = []

    last_reactions_ids = list(map(lambda x: x[1], previous_reactions_list))

    for reaction in reactions:

        if reaction.id not in last_reactions_ids:

            maximize = isMaximize(model, reaction, get_metabolite(model, metabolite_id), isReactant)

            if maximize is None:

                other_reactions_list.append((reaction, reaction.id, 'unknown',
                                             list_reactants_names(model, reaction.id),
                                             list_products_names(model, reaction.id)))

            elif maximize:

                reactions_list.append((reaction, reaction.id,
                                       simulate_reaction(model, reaction, maximize),
                                       list_reactants_names(model, reaction.id),
                                       list_products_names(model, reaction.id),
                                       get_reactants(model, reaction.id),
                                       get_products(model, reaction.id)))

            else:

                reactions_list.append((reaction, reaction.id,
                                       simulate_reaction(model, reaction, maximize),
                                       list_products_names(model, reaction.id),
                                       list_reactants_names(model, reaction.id),
                                       get_products(model, reaction.id),
                                       get_reactants(model, reaction.id)))

    return reactions_list, other_reactions_list

def get_reactions_by_role_fast(model, metabolite_id, isReactant, previous_reactions_list):

    reactions = get_reactions(model, metabolite_id)

    n_reactions = len(reactions)

    if n_reactions >= 20:

        reactions_list = []

        last_reactions_ids = list(map(lambda x: x[1], previous_reactions_list))

        n_test_set = round(n_reactions*0.1)

        reactions_in_list = list(reactions)

        test_set = reactions_in_list[:n_test_set]
        unknown_set = reactions_in_list[n_test_set:]

        for reaction in test_set:

            if reaction.id not in last_reactions_ids:

                maximize = isMaximize(model, reaction, get_metabolite(model, metabolite_id), isReactant)

                if maximize:

                    reactions_list.append((reaction, reaction.id,
                                           simulate_reaction(model, reaction, maximize),
                                           list_reactants_names(model, reaction.id),
                                           list_products_names(model, reaction.id),
                                           get_reactants(model, reaction.id),
                                           get_products(model, reaction.id)))

                else:

                    reactions_list.append((reaction, reaction.id,
                                           simulate_reaction(model, reaction, maximize),
                                           list_products_names(model, reaction.id),
                                           list_reactants_names(model, reaction.id),
                                           get_products(model, reaction.id),
                                           get_reactants(model, reaction.id)))

        for reaction in unknown_set:

            if not reaction.boundary and reaction.id not in last_reactions_ids:

                maximize = isMaximize(model, reaction, get_metabolite(model, metabolite_id), isReactant)

                if maximize:

                    reactions_list.append((reaction, reaction.id, 'unknown',
                                           list_reactants_names(model, reaction.id),
                                           list_products_names(model, reaction.id)))

                else:

                    reactions_list.append((reaction, reaction.id, 'unknown',
                                           list_products_names(model, reaction.id),
                                           list_reactants_names(model, reaction.id)))

        return reactions_list, []

    else:

        return get_reactions_by_role(model, metabolite_id, isReactant, previous_reactions_list)

def isMaximize(model, reaction, metabolite, isReactant):

    if float(reaction.upper_bound) > float(0.0) and float(reaction.lower_bound) < float(0.0):

        if isReactant:

            products = get_products(model, reaction.id)

            for product in products:

                if product.id == metabolite.id: return True

            return False

        else:

            reactants = get_reactants(model, reaction.id)

            for reactant in reactants:

                if reactant.id == metabolite.id: return True

            return False

    else:

        if float(reaction.upper_bound) > float(0.0):

            if isReactant:

                products = get_products(model, reaction.id)

                for product in products:

                    if product.id == metabolite.id: return True

                return None

            else:

                reactants = get_reactants(model, reaction.id)

                for reactant in reactants:

                    if reactant.id == metabolite.id: return True

                return None

        if float(reaction.lower_bound) < float(0.0):

            if isReactant:

                reactants = get_reactants(model, reaction.id)

                for reactant in reactants:

                    if reactant.id == metabolite.id: return False

                return None

            else:

                products = get_products(model, reaction.id)

                for product in products:

                    if product.id == metabolite.id: return False

                return None

@NodeCache
def simulate_reaction(model, reaction, isMaximize, tol=1E-06):

    with model as m:

        if isMaximize:

            set_objective_function(m, reaction.id)

            solution = m.slim_optimize()

            return evalSlimSol(solution, tol)

        else:

            get_reaction(m, reaction.id).bounds = (-999999, 0)

            set_objective_function(m, reaction.id)

            solution = m.optimize(objective_sense='minimize')

            return evalSol(solution, tol)

@NodeCache
def simulate_reactants(model, node, products, tol=1E-06):

    with model as m:

        for product in products:

            create_unbalenced_reaction(m, product.id, (-999999, 0))

        reaction_id = create_unbalenced_reaction(m, node.id, (0, 999999))

        set_objective_function(m, reaction_id)

        # using pFBA
        # solution = flux_analysis.pfba(m, reactions=[drain_name])

        # using fba, which is much much much faster than pFBA
        solution = m.slim_optimize()

    return evalSlimSol(solution, tol)

@NodeCache
def simulate_products(model, node, reactants, tol=1E-06):

    with model as m:

        for reactant in reactants:

            create_unbalenced_reaction(m, reactant.id, (0, 999999))

        reaction_id = create_unbalenced_reaction(m, node.id, (-999999, 0))

        set_objective_function(m, reaction_id)

        # using pFBA
        # solution = flux_analysis.pfba(m, reactions=[drain_name])

        # using fba, which is much much much faster than pFBA
        solution = m.optimize(objective_sense='minimize')

    return evalSol(solution, tol)

def evalSol(solution, tol):

    if np.isnan(solution.objective_value): return False

    if abs(solution.objective_value) < tol: return False

    if solution.status != 'optimal': return False

    return True

def evalSlimSol(solution, tol):

    if np.isnan(solution): return False

    if abs(solution) < tol: return False

    return True

def singleReactionKO(model, reaction_id, objective, exchange_prefix = None, tol = 1E-06):

    set_objective_function(model, reaction_id)

    initial_sol = model.optimize(objective_sense=objective)

    if not evalSol(initial_sol, tol):

        raise Exception("Objective value is zero. Model must have a valid solution objective value")

    # exchange reactions to remove
    if exchange_prefix is not None:

        reactions_remove = [i.id for i in model.reactions if exchange_prefix in i.name]

    else:

        reactions_remove = [i.id for i in model.reactions if 'Drainfor' in i.name or 'EX_' in i.id]

    for i in model.exchanges:

        if i.id not in reactions_remove:
            reactions_remove.append(i.id)

    reactions_remove.append(reaction_id)

    with model as m:
        reac_list = [get_reaction(m, i.id) for i in m.reactions if i.id not in reactions_remove]

    # multiprocessing is somehow taking more than just using a single core
    solution = single_reaction_deletion(model, reaction_list=reac_list, processes=1)

    solution.index = map(lambda x: x.__str__().replace('frozenset({\'', '').replace('\'})', ''), solution.index)

    solution.loc[:, 'growth'][np.isnan(solution.loc[:, 'growth'])] = float(0)
    solution.loc[:, 'growth'] = solution.loc[:, 'growth'].apply(lambda x: float(abs(x)))

    sr_lethal = solution[solution.loc[:, 'growth'] <= tol]

    if sr_lethal.shape[0] == 0:

        print(reaction_id, " failed", " since there is no KO available")
        print(reaction_id, " trying exchange reactions")

        reactions_remove = []

        reactions_remove.append(reaction_id)

        with model as m:
            reac_list = [get_reaction(m, i.id) for i in m.reactions if i.id not in reactions_remove]

        solution = single_reaction_deletion(model, reaction_list=reac_list, processes=1)

        solution.index = map(lambda x: x.__str__().replace('frozenset({\'', '').replace('\'})', ''), solution.index)

        solution.loc[:, 'growth'][np.isnan(solution.loc[:, 'growth'])] = float(0)
        solution.loc[:, 'growth'] = solution.loc[:, 'growth'].apply(lambda x: float(abs(x)))

        sr_lethal = solution[solution.loc[:, 'growth'] <= tol]

    if sr_lethal.shape[0] == 0:

        print(reaction_id, " failed", " since there is no KO available")

        return model

    reactions_to_select = list(sr_lethal.index)
    selected = random.choice(reactions_to_select)

    model.reactions.get_by_id(selected).bounds = (0.0, 0.0)

    set_objective_function(model, reaction_id)

    last_sol = model.optimize(objective_sense=objective)

    if evalSol(last_sol, tol):
        print(reaction_id, " failed", " since last solution is ", last_sol)
        raise Exception("The KO is not lethal, so cobrapy single reaction deletion solution is somehow incorrect")

    return model