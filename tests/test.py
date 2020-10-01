import os
import time
from unittest import TestCase, TestLoader, TextTestRunner

from tests.validation import biomass_model_processing
from bioiso import BioISO
from bioiso.wrappers.cobraWrapper import load, set_solver


class TestBioISO(TestCase):

    def setUp(self):
        self.startTime = time.time()

        self.model_name = 'iDS372'
        model_path = os.getcwd() + '/models/' + self.model_name + '.xml'
        self.reaction_to_eval = 'Biomass_assembly_C3_cytop'
        self.objective = 'maximize'
        solver = 'cplex'
        self.level = 2
        self.fast = False

        model = load(model_path)

        set_solver(model, solver)

        growth, self.m, reactions, metabolites = biomass_model_processing[self.model_name](model)

        self.presults = os.getcwd() + '/test_results/'

        if not os.path.exists(self.presults):
            os.mkdir(self.presults)

        with open(self.presults + 'GrowthCompoundRates' + self.model_name + '.txt', "w") as file:
            file.writelines('Growth: ' + str(growth) + '\n')
            file.writelines('Reactions: ' + str(reactions) + '\n')
            file.writelines('Metabolites: ' + str(metabolites) + '\n')

    def tearDown(self):
        t = time.time() - self.startTime
        print('%s: %.9f' % (self.id(), t))

    def test_BioISO(self):
        self.startTime = time.time()

        # # one Bioiso instance for each evaluation
        bio = BioISO(self.reaction_to_eval, self.m, self.objective)
        bio.run(self.level, self.fast)

        # # bio.write_results performs the get_tree and writes the results
        bio.write_results(self.presults + 'BioISOResults' + self.model_name + self.reaction_to_eval + '.json')

        assert bio.results['M_root_M_root_M_root_product']['analysis']


if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestBioISO)
    TextTestRunner(verbosity=0).run(suite)
