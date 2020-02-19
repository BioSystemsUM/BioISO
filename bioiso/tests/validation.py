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

    get_reaction(model, 'R00199_C3_cytop').bounds = (0.0, 0.0)

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

    citrate_rate = model.optimize().objective_value

    return ('EX_C00011_extr', citrate_rate), model

biomass_model_processing = {'iDS372' : iDS372,
                    'iJO1366' : iJO1366,
                    'iBsu1103' : iBsu1103,
                    'iTO977' : iTO977,
                    'iOD907' : iOD907}

compound_model_processing = {'iDS372' : iDS372_compound,
                    'iJO1366' : iJO1366_compound,
                    'iBsu1103' : iBsu1103_compound,
                    'iTO977' : iTO977_compound,
                    'iOD907' : iOD907_compound}

def read_and_processing_models(path, solver, biomass = True):

    print("Reading and processing models")

    models = {}

    if path:
        if not path.endswith('/'):
            path = path + '/'

    else:
        path = os.getcwd() + '/models/'

    for model_name in os.listdir(path):

        model = load(path + model_name)

        set_solver(model, solver)

        p = os.getcwd() + '/results/GrowthCompoundRates' + model_name[:-4] + '.txt'

        if biomass:

            rate, m = biomass_model_processing[model_name[:-4]](model)

            with open(p, "w") as file:
                file.writelines('Growth Rate: ' + str(rate) + '\n')

        else:

            rate, m = compound_model_processing[model_name[:-4]](model)

            with open(p, "w") as file:
                file.writelines('Compound Rate: ' + str(rate) + '\n')

        models[model_name[:-4]] = m

    print("Models are ready")

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

    for modelKey in models:

        print("Starting BioISO with model {}, reaction {} and objective {}".format(models[modelKey],
                                                                                   reactions[modelKey],
                                                                                   objectives[modelKey]))

        try:

            t0 = time.time()

            bio = Bioiso(reactions[modelKey], models[modelKey], objectives[modelKey])
            bio.run(level, fast)

            res = bio.get_tree()

            results[modelKey] = res

            path = os.getcwd() + '/results/'
            bio.write_results(path + fname + modelKey + reactions[modelKey] + objectives[modelKey] + '.json')

            t1 = time.time()

            print("BioISO has finished with running time of {}".format(str(t1-t0)))
            print("")


        except Exception as e:
            print("{} model failed".format(models[modelKey]))
            # raise e

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

def apply_ko(models, reactions, objectives):

    print("Applying KO")

    if not isinstance(models, dict):
        raise TypeError("models arg must be an {}".format(dict.__name__))

    if not isinstance(reactions, dict):
        raise TypeError("reactions arg must be an {}".format(dict.__name__))

    if not isinstance(objectives, dict):
        raise TypeError("objectives arg must be an {}".format(dict.__name__))

    if len(models) != len(reactions) or len(models) != len(objectives) or len(reactions) != len(objectives):
        raise ValueError("models, reactions and objectives args must have the same lenght. Current "
                         "dimensions models - {} reactions - {} objectives - {}".format(str(len(models)),
                                                                                        str(len(reactions)),
                                                                                        str(len(objectives))))

    new_models = {}

    for modelKey in models:

        new_models[modelKey] = singleReactionKO(models[modelKey], reactions[modelKey], objectives[modelKey])

    print("New models are ready")

    return new_models

def full_pipeline(biomass_reactions, biomass_objectives, biomass_fname, reactions_compounds, objectives_compounds, compounds_fname, solver = 'cplex', level =2 , fast = False):

    results, counts, results_compounds, counts_compounds = None, None, None, None

    if biomass_reactions:
        models = read_and_processing_models(False, solver)
        new_models = apply_ko(models, biomass_reactions, biomass_objectives)
        results = runBioiso_and_write_results(new_models, biomass_reactions, biomass_objectives, level, fast, biomass_fname + 'Results')
        counts = searchspacecounts(results, biomass_fname + 'Counts')

    if reactions_compounds:
        models = read_and_processing_models(False, solver, biomass=False)
        new_models_compounds = apply_ko(models, reactions_compounds, objectives_compounds)
        results_compounds = runBioiso_and_write_results(new_models_compounds, reactions_compounds, objectives_compounds, level, fast, compounds_fname + 'Results')
        counts_compounds = searchspacecounts(results_compounds, compounds_fname + 'Counts')

    return results, counts, results_compounds, counts_compounds

if __name__ == "__main__":

    solver = 'cplex'
    level = 2
    fast = False
    biomass_fname = 'BioISOBiomass'
    compounds_fname = 'BioISOCompounds'

    biomass_reactions = {
                 'iDS372': 'Biomass_assembly_C3_cytop',
                 'iJO1366' : 'Ec_biomass_iJO1366_core_53p95M',
                 'iOD907' : 'Biomass_cyto',
                 'iTO977' : 'CBIOMASS',
                 'iBsu1103' : 'bio00006'
                 }

    biomass_objectives = {
                  'iDS372': 'maximize',
                  'iJO1366' : 'maximize',
                  'iOD907' : 'maximize',
                  'iTO977' : 'maximize',
                  'iBsu1103' : 'maximize'
                  }

    reactions_compounds = {
                 'iDS372': 'R00703_C3_cytop',
                 'iJO1366' : 'ACKr',
                 'iOD907' : 'R00209_C4_mito',
                 'iTO977' : 'PDA1_2',
                 'iBsu1103' : 'rxn00227'
                 }

    objectives_compounds = {
                  'iDS372': 'maximize',
                  'iJO1366' : 'minimize',
                  'iOD907' : 'maximize',
                  'iTO977' : 'maximize',
                  'iBsu1103' : 'maximize'
                  }

    full_pipeline(biomass_reactions, biomass_objectives, biomass_fname, reactions_compounds, objectives_compounds, compounds_fname, solver, level, fast)



