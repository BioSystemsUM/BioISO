from bioiso.wrappers.cobraWrapper import load, set_solver, set_objective_function, get_reaction, singleReactionKO
from bioiso.core.bioiso import Bioiso
from bioiso.utils.bioisoUtils import searchSpaceSize, bioisoSearchSpace
import os
import json
import time
import warnings
warnings.filterwarnings("ignore")

def iDS372(model):

    set_objective_function(model, 'Biomass_assembly_C3_cytop')

    growth = model.optimize().objective_value

    set_objective_function(model, 'EX_C00186_extr')

    lactate_rate = model.optimize().objective_value

    return ('Biomass_assembly_C3_cytop', growth), ('EX_C00186_extr', lactate_rate), model

def iJO1366(model):

    set_objective_function(model, 'Ec_biomass_iJO1366_core_53p95M')

    get_reaction(model, 'EX_o2_LPAREN_e_RPAREN_').bounds = (0.0, 0.0)

    growth = model.optimize().objective_value

    set_objective_function(model, 'EX_ac_LPAREN_e_RPAREN_')

    acetate_rate = model.optimize().objective_value

    return ('Ec_biomass_iJO1366_core_53p95M', growth), ('EX_ac_LPAREN_e_RPAREN_', acetate_rate), model

def iBsu1103(model):

    set_objective_function(model, 'bio00006')

    for exchange in model.exchanges:
        exchange.lower_bound = -10.0

    growth = model.optimize().objective_value

    set_objective_function(model, 'rxn00225')

    acetate_rate = model.optimize().objective_value

    return ('bio00006', growth), ('rxn00225', acetate_rate), model

def iTO977(model):

    set_objective_function(model, 'CBIOMASS')

    model.reactions.get_by_id('GLCxtI').bounds = (-10.0, 10.0)

    growth = model.optimize().objective_value

    set_objective_function(model, 'EX_m11')

    co2_rate = model.optimize().objective_value

    return ('CBIOMASS', growth), ('EX_m11', co2_rate), model

def iOD907(model):

    set_objective_function(model, 'Biomass_cyto')

    model.reactions.get_by_id('EX_C00267_extr').lower_bound = -10.0

    growth = model.optimize().objective_value

    set_objective_function(model, 'EX_C00011_extr_b')

    co2_rate = model.optimize().objective_value

    return ('Biomass_cyto', growth), ('EX_C00011_extr_b', co2_rate), model

model_processing = {'iDS372' : iDS372,
                    'iJO1366' : iJO1366,
                    'iBsu1103' : iBsu1103,
                    'iTO977' : iTO977,
                    'iOD907' : iOD907}

def read_and_processing_models(solver, path = None):

    models = {}

    if not path:
        path = os.getcwd() + '/models/'
    else:
        if not path.endswith('/'):
            path = path+'/'

    for model_name in os.listdir(path):

        model = load(path + model_name)

        set_solver(model, solver)

        growth, compound, m = model_processing[model_name[:-4]](model)

        p = os.getcwd() + '/results/GrowthCompoundRates' + model_name[:-4] + '.txt'

        with open(p, "w") as file:
            file.writelines('Growth: ' + str(growth) + '\n')
            file.writelines('Compound: ' + str(compound))

        models[model_name[:-4]] = m

    return models

def runBioiso_and_write_results(models, reactions, objectives, level, fast, fname):

    if not isinstance(models, dict):
        raise TypeError("models arg must be an {}".format(dict.__name__))

    if not isinstance(reactions, dict):
        raise TypeError("reactions arg must be an {}".format(dict.__name__))

    if not isinstance(objectives, dict):
        raise TypeError("objectives arg must be an {}".format(dict.__name__))

    if not isinstance(level, int):
        raise TypeError("level arg must be an {}".format(int.__name__))

    if len(models) != len(reactions) or len(models) != len(objectives) or len(reactions) != len(objectives):
        raise ValueError("models, reactions and objectives args must have the same lenght. Current "
              "dimensions models - {} reactions - {} objectives - {}".format(str(len(models)), str(len(reactions)), str(len(objectives))))

    results = {}

    for modelKey, reactionKey, objectiveKey in zip(models, reactions, objectives):

        try:

            t0 = time.time()

            bio = Bioiso(reactions[reactionKey], models[modelKey], objectives[objectiveKey])
            bio.run(level, fast)

            res = bio.get_tree()

            results[modelKey] = res

            path = os.getcwd() + '/results/'
            bio.write_results(path + fname + modelKey + reactions[reactionKey] + objectives[objectiveKey] + '.json')

            t1 = time.time()
            print(t1-t0)

        except Exception as e:
            print("{} model failed".format(models[modelKey]))
            raise e

    return results

def load_results(path):

    with open(path) as jsonfile:
        return json.load(jsonfile)

def searchspacecounts(results, fname):

    res = {}

    for key in results:

        searchspacesize = searchSpaceSize(results[key])
        bioisosearchspacesize = bioisoSearchSpace(results[key])

        res[key + 'searchspacesize'] = searchspacesize
        res[key + 'bioisosearchspacesize'] = bioisosearchspacesize

        p = os.getcwd() + '/results/'+ fname + key + '.txt'

        with open(p, "w") as file:
            file.writelines('Total: ' + str(searchspacesize[0]) + '\n')
            file.writelines('Total reactions: ' + str(searchspacesize[1]) + '\n')
            file.writelines('Total metabolites: ' + str(searchspacesize[2]) + '\n')
            file.writelines('Unique total: ' + str(searchspacesize[3]) + '\n')
            file.writelines('Unique reactions: ' + str(searchspacesize[4]) + '\n')
            file.writelines('Unique metabolites: ' + str(searchspacesize[5]) + '\n')
            file.writelines('BioISO Total: ' + str(bioisosearchspacesize[0]) + '\n')
            file.writelines('BioISO Total reactions: ' + str(bioisosearchspacesize[1]) + '\n')
            file.writelines('BioISO Total metabolites: ' + str(bioisosearchspacesize[2]) + '\n')
            file.writelines('BioISO Unique total: ' + str(bioisosearchspacesize[3]) + '\n')
            file.writelines('BioISO Unique reactions: ' + str(bioisosearchspacesize[4]) + '\n')
            file.writelines('BioISO Unique metabolites: ' + str(bioisosearchspacesize[5]))

    return res

def apply_ko(models, reactions):

    if not isinstance(models, dict):
        raise TypeError("models arg must be an {}".format(dict.__name__))

    if not isinstance(reactions, dict):
        raise TypeError("reactions arg must be an {}".format(dict.__name__))

    if not isinstance(level, int):
        raise TypeError("level arg must be an {}".format(int.__name__))

    if len(models) != len(reactions):
        raise ValueError("models and reactions args must have the same lenght. Current "
              "dimensions models - {} reactions - {}".format(str(len(models)), str(len(reactions))))

    new_models = {}

    for modelKey, reactionsKey in zip(models, reactions):

        new_models[modelKey] = singleReactionKO(models[modelKey], reactions[reactionsKey])

    return new_models


if __name__ == "__main__":

    solver = 'cplex'
    level = 2
    fast = False

    # reactions = {
    #              'iDS372': 'Biomass_assembly_C3_cytop',
    #              'iJO1366' : 'Ec_biomass_iJO1366_core_53p95M',
    #              'iOD907' : 'Biomass_cyto',
    #              'iTO977' : 'CBIOMASS',
    #              'iBsu1103' : 'bio00006'
    #              }
    #
    # objectives = {
    #               'iDS372': 'maximize',
    #               'iJO1366' : 'maximize',
    #               'iOD907' : 'maximize',
    #               'iTO977' : 'maximize',
    #               'iBsu1103' : 'maximize'
    #               }

    # reactions_compounds = {
    #              'iDS372': 'EX_C00186_extr',
    #              'iJO1366' : 'EX_ac_LPAREN_e_RPAREN_',
    #              'iOD907' : 'EX_C00011_extr_b',
    #              'iTO977' : 'EX_m11',
    #              'iBsu1103' : 'rxn00225'
    #              }

    # objectives_compounds = {
    #               'iDS372': 'maximize',
    #               'iJO1366' : 'maximize',
    #               'iOD907' : 'maximize',
    #               'iTO977' : 'maximize',
    #               'iBsu1103' : 'maximize'
    #               }

    reactions = {
                 'iJO1366' : 'Ec_biomass_iJO1366_core_53p95M'
                 }

    objectives = {
                  'iDS372': 'maximize',
                  }

    models = read_and_processing_models(solver)

    new_models = apply_ko(models, reactions)

    results = runBioiso_and_write_results(new_models, reactions, objectives, level, fast, 'BioISOResults')

    counts = searchspacecounts(results, 'BioISOCounts')

