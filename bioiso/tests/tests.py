import os
import cobra
import numpy as np
import random
from cobra.test import create_test_model
from cobra.flux_analysis import single_reaction_deletion
from bioiso.wrappers.cobraWrapper import load

path = os.getcwd() + '/models/' + 'iDS372.xml'

m = load(path)

print(m.optimize())


reac_ids_to_remove = [i.id for i in m.reactions if 'Drainfor' in i.name or 'EX_' in i.id]

with m as modell:
    reac_list = [modell.reactions.get_by_id(i.id) for i in modell.reactions if i.id not in reac_ids_to_remove]

print(len(reac_list))
print(len(m.reactions))

sr = single_reaction_deletion(m, reaction_list = reac_list, processes=1)

sr.loc[:,'growth'][np.isnan(sr.loc[:,'growth'])] = float(0)
sr.loc[:, 'growth'] = sr.loc[:, 'growth'].apply(lambda x: float(abs(x)))
#cobra is weird as fuck
sr.index = map(lambda x: x.__str__().replace('frozenset({\'', '').replace('\'})', ''), sr.index)
tol = 1E-6
sr_letal = sr[sr.loc[:, 'growth'] <= 1E-6]
reactions_to_select = list(sr_letal.index)
selected = random.choice(reactions_to_select)

m.reactions.get_by_id(selected).bounds = (0.0, 0.0)

print(m.optimize())

# path = os.getcwd() + '/models/' + 'iDS372.xml'
# EX_C00186_extr lactate
# R_R00703_C3_cytop lactate dehydrogenase
# no o2 already
# glucose and wild-type conditions

# path = os.getcwd() + '/models/' + 'iJO1366.xml'
# EX_ac_LPAREN_e_RPAREN_ acetate
# R_ACKr_cytop acetate kinase
# no o2 EX_o2_LPAREN_e_RPAREN_
# -10 glucose and wild-type conditions

# path = os.getcwd() + '/models/' + 'iBsu1103.xml'
# EX_cpd00029_e acetate
# rxn00225 acetate kinase
# -10 all exchange reactions conditions

# path = os.getcwd() + '/models/' + 'iTO977.xml'
# EX_m11 CO2
# PDA1 Pyruvate dehydrogenase CO2
# KGD1 Alpha-ketoglutarate dehydrogenase CO2
# -10 m.reactions.get_by_id('GLCxtI').bounds glucose uptake with o2

# path = os.getcwd() + '/models/' + 'iOD907.xml'
# EX_C00011_extr_b CO2
# R00754_C3_cyto Acetaldehyde dehydrogenase Ethanol
# -10 m.reactions.get_by_id('EX_C00267_extr').lower_bound glucose uptake with o2



