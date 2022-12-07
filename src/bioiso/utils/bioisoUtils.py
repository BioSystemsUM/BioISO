import cobra
import time
from threading import Thread
import functools


class Node:

    def __init__(self, identifier=None, name=None, compartment=None, is_reactant=True):
        self.id = identifier
        self.name = name
        self.compartment = compartment
        self.is_reactant = is_reactant
        self.reactions_list = []
        self.other_reactions_list = []
        self.next = []
        self.previous = []
        self.analysis = None
        self.isLeaf = False

    def get_hash(self, stringify=False):

        role = 'product'
        if self.is_reactant:
            role = 'reactant'

        if stringify:
            return '_'.join((str(self.id), str(self.name), str(self.compartment), role))

        return self.id, self.name, self.compartment, role

    def get_next(self):
        return self.next

    def get_previous(self):
        return self.previous

    def has_next_nodes(self):
        return len(self.get_next()) > 0

    def has_next_by_hash(self, hash_tuple):

        next_nodes = self.get_next()

        for node in next_nodes:

            if node.get_hash() == hash_tuple:
                return node

        return

    def has_next(self, identifier):

        next_nodes = self.get_next()

        for node in next_nodes:

            if node.id == identifier:
                return node

        return

    def has_next_by_name(self, name):

        next_nodes = self.get_next()

        for node in next_nodes:

            if node.name == name:
                return node

        return


class NodeCache:
    bioiso_instances = {}

    @classmethod
    def new_node_cache(cls, instance):
        cls.bioiso_instances[instance] = {}
        return cls.bioiso_instances

    def __init__(self, func):
        self.function = func
        self._name = func.__name__

    def __call__(self, *args, **kwargs):

        bioiso_instance = args[0]

        if bioiso_instance not in self.bioiso_instances:
            analysis = self.function(*args, **kwargs)
            return analysis

        else:

            composed_id = self.create_composed_ids(self._name, args[1:])

            node_registry = self.bioiso_instances[bioiso_instance]

            if composed_id in node_registry:
                return node_registry[composed_id]

            else:
                analysis = self.function(*args, **kwargs)
                node_registry[composed_id] = analysis
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
                    composed_id = composed_id + str(met.id)

        return composed_id


def evaluate_side(boolean):
    if boolean:
        return 'Reactant'
    else:
        return 'Product'


def searchSpaceSizeRecursive(metabolites, reactions, tree):

    for key, values in tree.items():

        metabolites.add(values['identifier'])

        for reaction in values['reactions']:
            reactions.add(reaction[0])

        searchSpaceSizeRecursive(metabolites, reactions, values['next'])


def searchSpaceSize(tree):
    # length of all next nodes and reactions but without repeat

    metabolites = set()
    reactions = set()

    searchSpaceSizeRecursive(metabolites, reactions, tree)

    return len(reactions)+1, len(metabolites)-1


def bioisosearchSpaceSizeRecursive(metabolites, reactions, tree):

    for key, values in tree.items():

        if not values['analysis']:

            metabolites.add(values['identifier'])

        for reaction in values['reactions']:

            if not reaction[1]:
                reactions.add(reaction[0])

        bioisosearchSpaceSizeRecursive(metabolites, reactions, values['next'])


def bioisosearchSpaceSize(tree):
    # length of all next nodes and reactions but without repeat

    metabolites = set()
    reactions = set()

    bioisosearchSpaceSizeRecursive(metabolites, reactions, tree)

    return len(reactions), len(metabolites)


def timeit(function):
    def timed(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()

        print("### {} took {} ####".format(function.__name__, str(t1 - t0)))

        return result

    return timed


def timeout(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        seconds_before_timeout = args[0].timeout

        if seconds_before_timeout:

            res = [Exception('Bioiso timeout {} exceeded! Use Bioiso(timout=None) to bypass timeout'.format(
                seconds_before_timeout))]

            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as exc:
                    res[0] = exc

            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(seconds_before_timeout)
            except Exception as exc:
                print('error starting thread')
                raise exc
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret

            return ret

        else:

            return func(*args, **kwargs)

    return wrapper
