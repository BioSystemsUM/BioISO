import cobra
import time
from threading import Thread
import functools

class Node:

    def __init__(self, metabolite_id = None, name = None, is_reactant = True):
        self.id = metabolite_id
        self.name = name
        self.reactions_list = []
        self.other_reactions_list = []
        self.next = []
        self.previous = []
        self.flux = None
        self.isLeaf = False
        self.is_reactant = is_reactant

    def get_next(self):
        return self.next

    def get_previous(self):
        return self.previous

    def has_next_nodes(self):
        return self.get_next() == []

    def has_next(self, id):

        has = False

        next_nodes = self.get_next()

        for node in next_nodes:

            if node.id == id:
                return True, node

        return has, None

    def has_next_by_name(self, name):

        has = False

        next_nodes = self.get_next()

        for node in next_nodes:

            if node.name == name:
                return True, node

        return has, None

class NodeCache:

    node_registry = {}

    def __init__(self, func):
        self.function = func
        self._name = func.__name__

    def __call__(self, *args, **kwargs):

        composed_id = NodeCache.create_composed_ids(self._name, args)

        if composed_id in NodeCache.node_registry:
            return NodeCache.node_registry[composed_id]

        else:
            analysis = self.function(*args, **kwargs)
            NodeCache.node_registry[composed_id] = analysis
            return analysis

    @staticmethod
    def create_composed_ids(name, args):

        composed_id = '' + str(name)

        for arg in args:

            if isinstance(arg, cobra.core.model.Model):
                composed_id = composed_id + str(arg.name)

            if isinstance(arg, cobra.core.reaction.Reaction):
                composed_id = composed_id + str(arg.id)

            if isinstance(arg, Node):
                composed_id = composed_id + str(arg.id)

            if isinstance(arg, bool):
                composed_id = composed_id + str(arg)

            if isinstance(arg, list):

                for met in arg:

                    composed_id = composed_id +  str(met.id)

        return composed_id

def evaluate_side(bool):
    if bool:
        return 'Reactant'
    else:
        return 'Product'

def searchSpaceSize(results):

    # lenght of all next nodes and reactions but without repeat

    global metabolitesCache, reactionsCache, totalCache

    totalCache = {reaction[0]: 1 for reaction in results['M_fictitious']['reactions']}
    reactionsCache = {reaction[0]: 1 for reaction in results['M_fictitious']['reactions']}
    metabolitesCache = {}

    if len(totalCache) == 0:

        return 0, 0, 0, len(totalCache), len(reactionsCache), len(metabolitesCache), \
           totalCache, reactionsCache, metabolitesCache

    def searchSpaceSizeRecursive(results):

        global metabolitesCache, reactionsCache, totalCache

        for key in results:

            if key in metabolitesCache:
                metabolitesCache[key] += 1
                totalCache[key] += 1
            else:
                metabolitesCache[key] = 1
                totalCache[key] = 1

            for reaction in results[key]['reactions']:
                if reaction[0] in reactionsCache:
                    reactionsCache[reaction[0]] += 1
                    totalCache[reaction[0]] += 1
                else:
                    reactionsCache[reaction[0]] = 1
                    totalCache[reaction[0]] = 1

            if len(results[key]['next']) != 0:
                searchSpaceSizeRecursive(results[key]['next'])

            else:
                continue

    new_results = results['M_fictitious']['next']

    for key in new_results:

        if key in metabolitesCache:
            metabolitesCache[key] += 1
            totalCache[key] += 1
        else:
            metabolitesCache[key] = 1
            totalCache[key] = 1

        for reaction in new_results[key]['reactions']:
            if reaction[0] in reactionsCache:
                reactionsCache[reaction[0]] += 1
                totalCache[reaction[0]] += 1
            else:
                reactionsCache[reaction[0]] = 1
                totalCache[reaction[0]] = 1

        if len(new_results[key]['next']) != 0:

            searchSpaceSizeRecursive(new_results[key]['next'])

    totalTotal = 0
    for i in totalCache.values():
        totalTotal += i

    reactionsTotal = 0
    for i in reactionsCache.values():
        reactionsTotal += i

    metabolitesTotal = 0
    for i in metabolitesCache.values():
        metabolitesTotal += i

    return totalTotal, reactionsTotal, metabolitesTotal, len(totalCache), len(reactionsCache), len(metabolitesCache), \
           totalCache, reactionsCache, metabolitesCache


def bioisoSearchSpace(results):
    # length of the inherent nodes to a metabolite (precursor or successor) having a false analysis
    # but without the repeats

    global metabolitesCache, reactionsCache, totalCache

    totalCache = {reaction[0]: 1 for reaction in results['M_fictitious']['reactions']}
    reactionsCache = {reaction[0]: 1 for reaction in results['M_fictitious']['reactions']}
    metabolitesCache = {}

    if len(totalCache) == 0:

        return 0, 0, 0, len(totalCache), len(reactionsCache), len(metabolitesCache), \
           totalCache, reactionsCache, metabolitesCache

    def searchSpaceSizeRecursive(results):

        global metabolitesCache, reactionsCache, totalCache

        for key in results:

            if not results[key]['analysis']:

                if key in metabolitesCache:
                    metabolitesCache[key] += 1
                    totalCache[key] += 1
                else:
                    metabolitesCache[key] = 1
                    totalCache[key] = 1

                for reaction in results[key]['reactions']:

                    if not reaction[1]:

                        if reaction[0] in reactionsCache:
                            reactionsCache[reaction[0]] += 1
                            totalCache[reaction[0]] += 1
                        else:
                            reactionsCache[reaction[0]] = 1
                            totalCache[reaction[0]] = 1

                if len(results[key]['next']) != 0:
                    searchSpaceSizeRecursive(results[key]['next'])

                else:
                    continue

    new_results = results['M_fictitious']['next']

    for key in new_results:

        if not new_results[key]['analysis']:

            if key in metabolitesCache:
                metabolitesCache[key] += 1
                totalCache[key] += 1
            else:
                metabolitesCache[key] = 1
                totalCache[key] = 1

            for reaction in new_results[key]['reactions']:

                if not reaction[1]:

                    if reaction[0] in reactionsCache:
                        reactionsCache[reaction[0]] += 1
                        totalCache[reaction[0]] += 1
                    else:
                        reactionsCache[reaction[0]] = 1
                        totalCache[reaction[0]] = 1

            if len(new_results[key]['next']) != 0:
                searchSpaceSizeRecursive(new_results[key]['next'])

    totalTotal = 0
    for i in totalCache.values():
        totalTotal += i

    reactionsTotal = 0
    for i in reactionsCache.values():
        reactionsTotal += i

    metabolitesTotal = 0
    for i in metabolitesCache.values():
        metabolitesTotal += i

    return totalTotal, reactionsTotal, metabolitesTotal, len(totalCache), len(reactionsCache), len(metabolitesCache), \
           totalCache, reactionsCache, metabolitesCache


def timeit(function):
    def timed(*args, **kwargs):

        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()

        print("### {} took {} ####".format(function.__name__, str(t1-t0)))

        return result

    return timed

def timeout(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        seconds_before_timeout = args[0].__timeout__

        if seconds_before_timeout:

            res = [Exception('Bioiso timeout {} exceeded! Use Bioiso(timout=None) to bypass timeout'.format(
                seconds_before_timeout))]

            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e

            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(seconds_before_timeout)
            except Exception as e:
                print('error starting thread')
                raise e
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret

            return ret

        else:

            return func(*args, **kwargs)

    return wrapper
