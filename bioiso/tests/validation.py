from bioiso.wrappers.cobraWrapper import load, set_solver, set_objective_function, get_reaction, singleReactionKO
from bioiso.core.bioiso import Bioiso
from bioiso.utils.bioisoUtils import searchSpaceSize, bioisoSearchSpace
import os
import json
import time
import pandas as pd

import warnings

warnings.filterwarnings("ignore")


def iDS372(model):
    set_objective_function(model, 'Biomass_assembly_C3_cytop')

    growth = model.optimize().objective_value

    return ('Biomass_assembly_C3_cytop', growth), model


def iJO1366(model):
    set_objective_function(model, 'Ec_biomass_iJO1366_core_53p95M')

    get_reaction(model, 'EX_o2_LPAREN_e_RPAREN_').bounds = (0.0, 0.0)

    growth = model.optimize().objective_value

    return ('Ec_biomass_iJO1366_core_53p95M', growth), model


def iBsu1103(model):
    set_objective_function(model, 'bio00006')

    for exchange in model.exchanges:
        exchange.lower_bound = -10.0

    growth = model.optimize().objective_value

    return ('bio00006', growth), model


def iTO977(model):
    set_objective_function(model, 'CBIOMASS')

    get_reaction(model, 'GLCxtI').bounds = (-10.0, 10.0)

    growth = model.optimize().objective_value

    return ('CBIOMASS', growth), model


def iOD907(model):
    set_objective_function(model, 'Biomass_cyto')

    get_reaction(model, 'EX_C00267_extr').lower_bound = -10.0

    growth = model.optimize().objective_value

    return ('Biomass_cyto', growth), model


def iDS372_compound(model):
    set_objective_function(model, 'EX_C00186_extr')

    lactate_rate = model.optimize().objective_value

    return ('EX_C00186_extr', lactate_rate), model


def iJO1366_compound(model):
    set_objective_function(model, 'EX_ac_LPAREN_e_RPAREN_')

    get_reaction(model, 'EX_o2_LPAREN_e_RPAREN_').bounds = (0.0, 0.0)

    acetate_rate = model.optimize().objective_value

    return ('EX_ac_LPAREN_e_RPAREN_', acetate_rate), model


def iBsu1103_compound(model):
    set_objective_function(model, 'EX_cpd00029_e')

    for exchange in model.exchanges:
        exchange.lower_bound = -10.0

    get_reaction(model, 'rxn00152').bounds = (0.0, 0.0)
    get_reaction(model, 'rxn00173').bounds = (0.0, 0.0)

    acetate_rate = model.optimize().objective_value

    return ('EX_cpd00029_e', acetate_rate), model


def iTO977_compound(model):
    set_objective_function(model, 'CO2xtO')

    get_reaction(model, 'GLCxtI').bounds = (-10.0, 10.0)

    co2_rate = model.optimize().objective_value

    return ('CO2xtO', co2_rate), model


def iOD907_compound(model):
    set_objective_function(model, 'EX_C00158_extr')

    get_reaction(model, 'EX_C00267_extr').lower_bound = -10.0

    get_reaction(model, 'R00355_C3_cyto').bounds = (0.0, 0.0)
    get_reaction(model, 'R00344_C3_cyto').bounds = (0.0, 0.0)
    get_reaction(model, 'T01268_C4_mito').bounds = (0.0, 0.0)
    get_reaction(model, 'R01731_C3_cyto').bounds = (0.0, 0.0)

    citrate_rate = model.optimize().objective_value

    return ('EX_C00158_extr', citrate_rate), model


biomass_model_processing = {'iDS372': iDS372,
                            'iJO1366': iJO1366,
                            'iBsu1103': iBsu1103,
                            'iTO977': iTO977,
                            'iOD907': iOD907}

compound_model_processing = {'iDS372': iDS372_compound,
                             'iJO1366': iJO1366_compound,
                             'iBsu1103': iBsu1103_compound,
                             'iTO977': iTO977_compound,
                             'iOD907': iOD907_compound}


def read_and_processing_models(path, solver, biomass=True):

    print("Reading and processing models")

    models = {}

    if path:
        if not path.endswith('/'):
            path = path + '/'

    else:
        path = os.getcwd() + '/models/'

    for model_name in os.listdir(path):

        if model_name == 'results':
            continue

        print("Reading {} model".format(model_name))

        model = load(path + model_name)

        set_solver(model, solver)

        p = os.getcwd() + '/results/GrowthCompoundRates' + model_name[:-4] + '.txt'

        if biomass:

            rate, m = biomass_model_processing[model_name[:-4]](model)

            with open(p, "a") as file:
                file.writelines('Growth Rate: ' + str(rate) + '\n')

        else:

            rate, m = compound_model_processing[model_name[:-4]](model)

            with open(p, "a") as file:
                file.writelines('Compound Rate: ' + str(rate) + '\n')

        models[model_name[:-4]] = m

    print("Models are ready")

    return models


def pipeline(models_path, reactions, objectives, biomass=True, solver='cplex', level=2, fast=False):

    if not isinstance(reactions, dict):
        raise TypeError("reactions arg must be an {}".format(dict.__name__))

    if not isinstance(objectives, dict):
        raise TypeError("objectives arg must be an {}".format(dict.__name__))

    if not isinstance(level, int):
        raise TypeError("level arg must be an {}".format(int.__name__))

    models = read_and_processing_models(models_path, solver, biomass=biomass)

    if len(models) != len(reactions) or len(models) != len(objectives) or len(reactions) != len(objectives):
        raise ValueError("models, reactions and objectives args must have the same lenght. "
                         "Current dimensions: "
                         "models - {} "
                         "reactions - {} "
                         "objectives - {}".format(str(len(models)), str(len(reactions)), str(len(objectives))))

    results_path = os.path.join(models_path, 'results')

    if not os.path.exists(results_path):
        os.mkdir(results_path)

    models_kos = {}

    print()
    print('Getting KOs')

    for modelKey in models:

        try:

            print("Getting KOs for {} model".format(modelKey))

            model_kos, model_all_kos = singleReactionKO(models[modelKey], reactions[modelKey], objectives[modelKey])

            models_kos[modelKey] = model_kos

            kos_df = pd.DataFrame(model_all_kos)

            f_name = results_path + '\{}_{}_{}_kos.csv'.format(modelKey, reactions[modelKey], objectives[modelKey])

            kos_df.to_csv(f_name)

        except ValueError as e:

            print("KO Error for {} model with exception {}".format(modelKey, e))

    print()
    print('KOs are ready for {} models'.format(str(len(models_kos))))

    trees = {}
    results = {}

    for modelKey in models_kos:

        trees[modelKey] = {}
        results[modelKey] = {}

        for ko in models_kos[modelKey]:

            print("Starting BioISO with model {}, "
                  "reaction {} "
                  "and objective {} "
                  "for ko {}".format(modelKey, reactions[modelKey], objectives[modelKey], ko))

            results[modelKey][ko] = {}

            t0 = time.time()

            with models[modelKey] as m:

                m.reactions.get_by_id(ko).bounds = (0.0, 0.0)

                bio = Bioiso(reactions[modelKey], m, objectives[modelKey])
                bio.run(level, fast)

                t1 = time.time()

                tree = bio.get_tree()

                trees[modelKey][ko] = tree

                f_name = results_path + '\{}_{}_{}_{}.json'.format(modelKey, ko, reactions[modelKey],
                                                                  objectives[modelKey])

                bio.write_results(f_name)

                results[modelKey][ko].update({'time': float(t1 - t0)})

            print("BioISO has finished with running time of {}".format(str(t1 - t0)))
            print("")

    print()
    print('Getting BioISO statistics for {} models'.format(str(len(trees))))

    for modelKey in trees:

        for ko in trees[modelKey]:

            searchspacesize = searchSpaceSize(trees[modelKey][ko])
            bioisosearchspacesize = bioisoSearchSpace(trees[modelKey][ko])

            results[modelKey][ko].update({
                'total': searchspacesize[0],
                'total_reactions': searchspacesize[1],
                'total_metabolites': searchspacesize[2],
                'unique': searchspacesize[3],
                'unique_reactions': searchspacesize[4],
                'unique_metabolites': searchspacesize[5],

                'bioiso_total': bioisosearchspacesize[0],
                'bioiso_total_reactions': bioisosearchspacesize[1],
                'bioiso_total_metabolites': bioisosearchspacesize[2],
                'bioiso_unique': bioisosearchspacesize[3],
                'bioiso_unique_reactions': bioisosearchspacesize[4],
                'bioiso_unique_metabolites': bioisosearchspacesize[5],
            })

    print()
    print('Writing ...')

    f_name = '_'.join([key[0:5] if len(key) > 5 else key for key in reactions.values()])

    f_name = results_path + '/' + f_name + '_bioiso.xlsx'

    dfs = []

    with pd.ExcelWriter(f_name) as writer:

        for modelKey in results:

            rows = list(results[modelKey].keys())
            cols = list(results[modelKey][rows[0]].keys())
            data = [[metric for metric in results[modelKey][ko].values()] for ko in results[modelKey]]

            df = pd.DataFrame(data=data, index=rows, columns=cols)
            df.to_excel(writer, sheet_name=modelKey)
            dfs.append(df)

    return dfs, results, trees


if __name__ == "__main__":

    models_path = os.getcwd()

    models_path = os.path.join(models_path, ('models'))

    solver = 'cplex'
    level = 2
    fast = False

    biomass_reactions = {
        'iDS372': 'Biomass_assembly_C3_cytop',
        'iJO1366': 'Ec_biomass_iJO1366_core_53p95M',
        'iOD907': 'Biomass_cyto',
        'iTO977': 'CBIOMASS',
        'iBsu1103': 'bio00006'
    }

    biomass_objectives = {
        'iDS372': 'maximize',
        'iJO1366': 'maximize',
        'iOD907': 'maximize',
        'iTO977': 'maximize',
        'iBsu1103': 'maximize'
    }

    compounds_reactions = {
        'iDS372': 'R00703_C3_cytop',
        'iJO1366': 'ACKr',
        'iOD907': 'R00351_C3_cyto',
        'iTO977': 'PDA1_2',
        'iBsu1103': 'rxn00227'
    }

    compounds_objectives = {
        'iDS372': 'maximize',
        'iJO1366': 'minimize',
        'iOD907': 'maximize',
        'iTO977': 'maximize',
        'iBsu1103': 'maximize'
    }

    bio_df = pipeline(models_path, biomass_reactions, biomass_objectives, biomass=True, solver='cplex', level=2, fast=False)

    compounds_df = pipeline(models_path, compounds_reactions, compounds_objectives, biomass=False, solver='cplex', level=2,
                            fast=False)
