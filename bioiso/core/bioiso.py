import json
from bioiso.wrappers.cobraWrapper import get_products, get_reactants, get_reaction, get_metabolite, \
    list_reactants_names, list_products_names, simulate_reaction, get_reactions_by_role, simulate_reactants, \
    simulate_products, get_reactions_by_role_fast
from bioiso.utils.bioisoUtils import Node, NodeCache, evaluate_side, timeout

class Bioiso:

    def __init__(self, reaction_id, model, objective_direction, fast = False, timeout=600):

        self.reaction_id = reaction_id
        self.model = model
        self.root = None
        self.results = None

        self.setObjective(objective_direction)

        # controlling cache
        # any time Bioiso starts the registry must be cleaned
        NodeCache.node_registry = {}

        # controlling recursion timeout
        self.__timeout__ = timeout

        # fast version ignores the evaluation of 90% of all reactions (by selecting the first ones) associated with a
        # node if this node has more than 20 reactions
        # the evaluation result of the remaining is unknown
        # also the node is set as leaf, and Bioiso doesn't go any further
        self.__fast__ = fast

        self.__principles_verified__ = False


    def changeTimeout(self, timeout):
        self.__timeout__ = timeout

    def setObjective(self, objective_direction):

        if objective_direction == 'maximize':
            self.__objective_direction__ = 'maximize'

        elif objective_direction == 'minimize':
            self.__objective_direction__ = 'minimize'

        else:
            try:
                print("Oops! {} is not a valid option! Please try maximize or minimize".format(str(objective_direction)))
                raise ValueError
            except TypeError as e:
                print("Please try maximize or minimize as an {} object".format(str(str.__name__)))
                raise e

    def __verify_reaction_principles__(self):

        '''Verifies the principles for testing the precursors of the input reaction
        Return Boolean (it has all the principles), products
        '''

        try:
            if self.__objective_direction__ == 'maximize':

                # get the products of the reaction
                products = get_products(self.model, self.reaction_id)

                self.__metabolite_is_reactant__ = False

                # creating the root metabolite
                self.__metabolite_id__ = 'M_fictitious'
                self.__metabolite_name__ = 'M_fictitious'

            else:

                # get the products of the reaction
                reactants = get_reactants(self.model, self.reaction_id)

                self.__metabolite_is_reactant__ = True

                # creating the root metabolite
                self.__metabolite_id__ = 'M_fictitious'
                self.__metabolite_name__ = 'M_fictitious'


            self.__principles_verified__ = True

        except KeyError as e:
            print("Please select a valid reaction identifier")
            raise e

        except:
            print("Unexpected error")
            raise

    def set_root(self):

        '''Setting the root, namely the product metabolite in the reaction provided
        If the model doesn't have a single product,
        a fictitious metabolite is set as root.
        It does not add this fictitious metabolite to the model though

        This is controlled in the __verify_reaction_principles__ method'''

        if not self.__principles_verified__: self.__verify_reaction_principles__()

        #setting the root
        self.root = Node(self.__metabolite_id__, self.__metabolite_name__, self.__metabolite_is_reactant__)

        reaction = get_reaction(self.model, self.reaction_id)

        # testing first the root reaction
        if self.__objective_direction__ == 'maximize':

            reaction_flux = simulate_reaction(self.model, reaction, True)

            self.root.reactions_list = [(reaction, reaction.id, reaction_flux,
                                         list_reactants_names(self.model, reaction.id),
                                         list_products_names(self.model, reaction.id),
                                        get_reactants(self.model, reaction.id),
                                        get_products(self.model, reaction.id))]

        else:

            reaction_flux = simulate_reaction(self.model, reaction, False)

            self.root.reactions_list = [(reaction, reaction.id, reaction_flux,
                                         list_products_names(self.model, reaction.id),
                                         list_reactants_names(self.model, reaction.id),
                                        get_products(self.model, reaction.id),
                                        get_reactants(self.model, reaction.id))]

        self.root.flux = reaction_flux

    def create_reactions_list(self, node, last_reaction_list, reactant):

        '''Creates the node reactions_list
        It populates the property reactions_list of the object node which is given as input
        with the reactions associated to the metabolite in the model
        the previous_reactions_list argument controls the appearances of repetitions.
        For instance, the metabolite protein should have only one reactions_list, the R_Protein.
        To do so, the R_Biomass should not be considered, even if the protein metabolite
        takes part into this reaction.'''

        # the reactions_list of this metabolite should not contain the previous reactions_list
        # Besides, if the metabolite is a reactant, only reactions where the metabolite takes part as product should be selected for reactions_list
        # Otherwise, reactions where the metabolite takes part as reactant should be selected for reactions_list

        if self.__fast__:
            node.reactions_list, node.other_reactions_list = get_reactions_by_role_fast(self.model, node.id, isReactant=
            reactant, previous_reactions_list=last_reaction_list)

        else:
            node.reactions_list, node.other_reactions_list = get_reactions_by_role(self.model, node.id, isReactant =
            reactant, previous_reactions_list=last_reaction_list)

    def create_next_nodes(self, node, leaf=False):

        if self.__fast__:

            if len(node.reactions_list) >= 20:

                node.isLeaf = True

            else:

                # let's iterate the property reactions_list of the node
                for reaction in node.reactions_list:
                    # for each reaction let's set the next nodes, namely the reactants of each reaction in reactions_list
                    self.__create_next_nodes__(node, reaction, leaf)

        else:

            # let's iterate the property reactions_list of the node
            for reaction in node.reactions_list:

                # for each reaction let's set the next nodes, namely the reactants of each reaction in reactions_list
                self.__create_next_nodes__(node, reaction, leaf)

    def __create_next_nodes__(self, node, reaction, leaf):

        reactants = reaction[5]
        products = reaction[6]

        for reactant in reactants:

            isNode, nextNode = node.has_next(reactant.id)

            if not isNode:

                #but if a metabolite with an equal name already exists,
                # it means that the metabolite is on a different compartment
                # let's create a composed name

                isNode2, _ = node.has_next_by_name(reactant.name)

                if isNode2:
                    composed_name = reactant.name + ' ' + reactant.compartment

                else:
                    composed_name = reactant.name


                #for each reactant in those reactions
                #create a new node using the id and name
                next_node = Node(reactant.id, composed_name, True)

                #the previous engine is similar to the next
                #the previous node of the next nodes is the current node itself
                next_node.previous = [node]

                #flag to control leafs
                #by default nodes are not leafs (see create_next_nodes)
                #populate tree controls the leafs
                next_node.isLeaf = leaf

                previous_reactions_list = node.reactions_list
                self.create_reactions_list(next_node, last_reaction_list = previous_reactions_list, reactant = True)

                #append next nodes of the current node, namely the reactants of the reactions where the node takes part into
                node.next.append(next_node)

                # the products in the reactions_list reaction are the ones that might impair the flux
                next_node.flux = simulate_reactants(self.model, next_node, products)

            else:

                if not nextNode.is_reactant:

                    composed_name = reactant.name + ' ' + reactant.compartment + ' as reactant'

                    next_node = Node(reactant.id, composed_name, True)
                    next_node.previous = [node]
                    next_node.isLeaf = leaf
                    previous_reactions_list = node.reactions_list
                    self.create_reactions_list(next_node, last_reaction_list=previous_reactions_list, reactant=True)
                    node.next.append(next_node)
                    next_node.flux = simulate_reactants(self.model, next_node, products)

        #products version
        for product in products:

            isNode, nextNode = node.has_next(product.id)

            if not isNode:

                isNode2, _ = node.has_next_by_name(product.name)

                if isNode2:
                    composed_name = product.name + ' ' + product.compartment

                else:
                    composed_name = product.name

                next_node = Node(product.id, composed_name, False)
                next_node.previous = [node]
                next_node.isLeaf = leaf

                previous_reactions_list = node.reactions_list
                self.create_reactions_list(next_node, last_reaction_list = previous_reactions_list, reactant = False)
                node.next.append(next_node)
                next_node.flux = simulate_products(self.model, next_node, reactants)

            else:

                if nextNode.is_reactant:

                    composed_name = product.name + ' ' + product.compartment + ' as product'

                    next_node = Node(product.id, composed_name, False)
                    next_node.previous = [node]
                    next_node.isLeaf = leaf
                    previous_reactions_list = node.reactions_list
                    self.create_reactions_list(next_node, last_reaction_list=previous_reactions_list, reactant=False)
                    node.next.append(next_node)
                    next_node.flux = simulate_products(self.model, next_node, reactants)

    def run(self, levels, fast = False):

        self.__fast__ = fast

        if not self.root:
            self.set_root()

        if self.results is not None:

            self.results = None

        return self.populate_tree(levels)

    @timeout
    def populate_tree(self, levels):

        '''
        create nodes associated with the root (biomass metabolite), namely each macromolecule or subunit,
        and subsequently
        '''

        self.levels = levels

        if self.levels < 0:
            return

        elif self.levels == 0:

            self.root.isLeaf = True

        elif self.levels == 1:

            self.create_next_nodes(self.root, leaf=True)

        else:
            # if level is higher than 2
            # a recursive function is called

            self.__populate_tree__([self.root], level=0)

    def __populate_tree__(self, nodes, level):

        '''Recursive hidden method to populate the tree after the level 1
        After level 1 this function is called with the next nodes of root, namely the first next nodes
        for each node in the list of next nodes,
        the recursive function goes deep until it reaches the maximum level (controlled by level == self.levels)
        going through each level in a given node,
        the recursive function creates the next_nodes of this node and then get these for next recursion
        in the end, the final nodes are set as leafs
        '''

        if level == self.levels:

            # exit for recursive function populate_sub_tree

            #each next node in the last level is set as leaf

            for node in nodes:
                node.isLeaf = True

            # finally, the recursion is ended
            return


        elif len(nodes) == 0:

            # exit for recursive function populate_sub_tree

            #each node without next node is set as leaf

            for node in nodes:
                node.isLeaf = True

            # finally, the recursion is ended
            return

        else:

            # if the recursive function has not reached the maximum level

            # iteration over the next nodes of the root and following next nodes of these
            for node in nodes:

                self.create_next_nodes(node)

                nodes = node.get_next()

                self.__populate_tree__(nodes, level+1)

    def get_tree(self):

        self.results = {self.root.name: {'analysis': self.root.flux,
                                        'role': evaluate_side(self.root.is_reactant),
                                        'reactions': [(child[1], child[2], child[3], child[4]) for child in
                                              self.root.reactions_list],
                                        'other_reactions': [(None, None, None, None)],
                                        'next': {}}}

        self.__get_tree__([self.root], self.results, level=0)

        return self.results

    def __get_tree__(self, nodes, dictionary, level):


        if level == self.levels+1:
            return

        else:

            for node in nodes:

                next_nodes = node.get_next()

                dictionary[node.name]['next'] = {
                    next_node.name: {'analysis': next_node.flux,
                                     'role': evaluate_side(next_node.is_reactant),
                                     'reactions': [(child[1], child[2], child[3], child[4]) for child in next_node.reactions_list],
                                     'other_reactions': [(child[1], child[2], child[3], child[4]) for child in
                                                   next_node.other_reactions_list],
                                     'next': {}} for next_node in next_nodes}

                self.__get_tree__(next_nodes, dictionary[node.name]['next'], level+1)

    def write_results(self, results_fname):

        if not self.results:

            self.get_tree()

        with open(results_fname, "w") as jsonfile:
            json.dump(self.results, jsonfile)