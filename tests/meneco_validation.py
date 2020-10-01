import os
import time

import pandas as pd
from meneco import run_meneco
from meneco.meneco import query, sbml


def read_and_processing_models(path, models_names):
    print("Reading and processing models")

    models = {}

    if path:
        if not path.endswith('/'):
            path = path + '/'

    else:
        path = os.getcwd() + '/meneco_models/'

    for model_name in os.listdir(path):

        if model_name in models_names:

            models[model_name] = {}

            print("Reading {} models".format(model_name))

            seeds = sbml.readSBMLseeds(os.path.join(path, model_name, model_name + '_seeds.xml'))
            targets = sbml.readSBMLtargets(os.path.join(path, model_name, model_name + '_targets.xml'))

            for model_ko in os.listdir(os.path.join(path, model_name)):

                if not model_ko.endswith('_seeds.xml') and not model_ko.endswith('_targets.xml'):

                    model = sbml.readSBMLnetwork(os.path.join(path, model_name, model_ko), 'draft')

                    models[model_name][model_ko] = (model, seeds, targets)

    print("Models are ready")

    return models


def pipeline(paths, models_names):

    if not isinstance(models_names, list):
        raise TypeError("models names arg must be an {}".format(list.__name__))

    results_path = os.path.join(os.getcwd(), 'meneco_validation_results')

    if not os.path.exists(results_path):
        os.mkdir(results_path)

    models = read_and_processing_models(paths, models_names)

    if len(models) != len(models_names):
        raise ValueError("models and names args must have the same length. "
                         "Current dimensions: "
                         "models - {} "
                         "models names - {} ".format(str(len(models)), str(len(models_names))))

    print()
    print('Getting Meneco statistics for {} models'.format(str(len(models))))

    results = {}

    for modelKey in models:

        results[modelKey] = {}

        for ko in models[modelKey]:

            draftnet, seeds, targets = models[modelKey][ko]
            t0 = time.time()
            model = query.get_unproducible(draftnet, targets, seeds)
            unproducible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
            t1 = time.time()

            results[modelKey][ko] = {
                'time': t1-t0,
                'total metabolites': len(targets),
                'meneco metabolites': len(unproducible),
            }

    print()
    print('Writing ...')

    f_name = '_'.join([key[0:5] if len(key) > 5 else key for key in models_names])

    f_name = results_path + '/' + f_name + '_meneco.xlsx'

    dfs = []

    with pd.ExcelWriter(f_name) as writer:

        for modelKey in results:
            rows = list(results[modelKey].keys())
            cols = list(results[modelKey][rows[0]].keys())
            data = [[metric for metric in results[modelKey][ko].values()] for ko in results[modelKey]]

            df = pd.DataFrame(data=data, index=rows, columns=cols)
            df.to_excel(writer, sheet_name=modelKey)
            dfs.append(df)

    print()
    print('Ready ...')

    return dfs, results


if __name__ == "__main__":

    import json
    import warnings

    warnings.filterwarnings("ignore")

    models_path = os.getcwd()

    models_path = os.path.join(models_path, 'meneco_models')

    # models_names = ['iDS372', 'iJO1366']

    # bio_df = pipeline(models_path, models_names)

    # filling

    # result = run_meneco(draftnet=os.path.join(models_path, 'iJO1366', 'SO4tex.xml'),
    #                     seeds=os.path.join(models_path, 'iJO1366', 'iJO1366_seeds.xml'),
    #                     targets=os.path.join(models_path, 'iJO1366', 'iJO1366_targets.xml'),
    #                     repairnet=os.path.join(models_path, 'iJO1366', 'bigg_universal_model.xml'),
    #                     enumeration=True,
    #                     json=True)
    #
    # with open('meneco_ijo1366_gap_fill_default.json', "w") as jsonfile:
    #     json.dump(result, jsonfile)

    result = run_meneco(draftnet=os.path.join(models_path, 'iDS372', 'R04568_C3_cytop.xml'),
                        seeds=os.path.join(models_path, 'iDS372', 'iDS372_seeds.xml'),
                        targets=os.path.join(models_path, 'iDS372', 'iDS372_targets.xml'),
                        repairnet=os.path.join(models_path, 'iDS372', 'kegg_universal_model.xml'),
                        enumeration=True,
                        json=True)

    with open('meneco_ids372_gap_fill_default.json', "w") as jsonfile:
        json.dump(result, jsonfile)

    result = run_meneco(draftnet=os.path.join(models_path, 'iJO1366', 'SO4tex.xml'),
                        seeds=os.path.join(models_path, 'iJO1366', 'iJO1366_seeds.xml'),
                        targets=os.path.join(models_path, 'iJO1366', 'iJO1366_bioiso_targets.xml'),
                        repairnet=os.path.join(models_path, 'iJO1366', 'bigg_universal_model.xml'),
                        enumeration=True,
                        json=True)

    with open('meneco_ijo1366_gap_fill_bioiso.json', "w") as jsonfile:
        json.dump(result, jsonfile)

    result = run_meneco(draftnet=os.path.join(models_path, 'iDS372', 'R04568_C3_cytop.xml'),
                        seeds=os.path.join(models_path, 'iDS372', 'iDS372_seeds.xml'),
                        targets=os.path.join(models_path, 'iDS372', 'iDS372_bioiso_targets.xml'),
                        repairnet=os.path.join(models_path, 'iDS372', 'kegg_universal_model.xml'),
                        enumeration=True,
                        json=True)

    with open('meneco_ids372_gap_fill_bioiso.json', "w") as jsonfile:
        json.dump(result, jsonfile)
