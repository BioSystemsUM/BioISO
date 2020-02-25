import os
from .validation import biomass_model_processing, compound_model_processing
from bioiso.core.bioiso import Bioiso
from bioiso.wrappers.cobraWrapper import load, set_solver

model_name = 'iDS372'
model_path = os.getcwd() + '/models/' + model_name + '.xml'
reaction_to_eval = 'Biomass_assembly_C3_cytop'
objective = 'maximize'
solver = 'cplex'
level = 2
fast = False

model = load(model_path)

set_solver(model, solver)

growth, m = biomass_model_processing[model_name](model)

presults = os.getcwd() + '/results/'

with open(presults + 'GrowthCompoundRates' + model_name + '.txt', "w") as file:
    file.writelines('Growth: ' + str(growth) + '\n')

#one Bioiso instance for each evaluation
bio = Bioiso(reaction_to_eval, model, objective)
bio.run(level, fast)

# getting the results
# res = bio.get_tree()

# bio.write_results is the get_tree and writes the results
bio.write_results(presults + 'BioISOResults' + model_name + reaction_to_eval + '.json')
