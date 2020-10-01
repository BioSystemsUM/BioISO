import json
from bioiso.wrappers.cobraWrapper import get_products, get_reactants, get_reaction, \
    list_reactants_ids, list_products_ids, simulate_reaction, get_reactions_by_role, simulate_reactants, \
    simulate_products, get_reactions_by_role_fast
from bioiso.utils.bioisoUtils import Node, NodeCache, evaluate_side, timeout


class BioISO:

    def __init__(self, reaction_id, model, objective_direction, fast=False, time_out=900):

        self.levels = 0

        self.reaction_id = reaction_id
        self.model = model
        self.root = None
        self.summary = {}
        self.results = {}

        self.objective_direction = 'maximize'
        self.setObjective(objective_direction)

        self.__verify_reaction_principles()

        # controlling cache
        # any time Bioiso starts the registry must be cleaned
        self.__id = self.reaction_id + '_' + self.model.id + '_' + self.objective_direction + '_' + str(id(self))
        self.nodes_cache = NodeCache.new_node_cache(self.__id)

        # controlling recursion timeout
        self.timeout = time_out

        # fast version ignores the evaluation of 90% of all reactions (by selecting the first ones) associated with a
        # node if this node has more than 20 reactions
        # the evaluation result of the remaining is unknown
        # also the node is set as leaf, and Bioiso doesn't go any further
        self.fast = fast

        self.__principles_verified = False

    def changeTimeout(self, timeout_time):
        self.timeout = timeout_time

    def setObjective(self, objective_direction):

        if objective_direction == 'maximize':
            self.objective_direction = 'maximize'

        elif objective_direction == 'minimize':
            self.objective_direction = 'minimize'

        else:
            try:
                print(
                    "Oops! {} is not a valid option! Please try maximize or minimize".format(str(objective_direction)))
                raise ValueError
            except TypeError as e:
                print("Please try maximize or minimize as an {} object".format(str(str.__name__)))
                raise e

    def __verify_reaction_principles(self):

        """Verifies the principles for testing the precursors of the input reaction
        Return Boolean (it has all the principles), products
        """

        try:
            if self.objective_direction == 'maximize':

                # get the products of the reaction
                products = get_products(self.model, self.reaction_id)

                self.__metabolite_is_reactant = False

                # creating the root metabolite
                self.__metabolite_id = 'M_root'
                self.__metabolite_name = 'M_root'

            else:

                # get the products of the reaction
                reactants = get_reactants(self.model, self.reaction_id)

                self.__metabolite_is_reactant = True

                # creating the root metabolite
                self.__metabolite_id = 'M_root'
                self.__metabolite_name = 'M_root'

            self.__principles_verified = True

        except KeyError as e:
            print("Please select a valid reaction identifier")
            raise e

        except BaseException:
            print("Unexpected error")
            raise

    def set_root(self):

        """Setting the root, namely the product metabolite in the reaction provided
        If the model doesn't have a single product,
        a fictitious metabolite is set as root.
        It does not add this fictitious metabolite to the model though

        This is controlled in the __verify_reaction_principles__ method"""

        if not self.__principles_verified:
            self.__verify_reaction_principles()

        # setting the root
        self.root = Node(identifier=self.__metabolite_id,
                         name=self.__metabolite_name,
                         compartment=self.__metabolite_name,
                         is_reactant=self.__metabolite_is_reactant)

        reaction = get_reaction(self.model, self.reaction_id)

        # testing first the root reaction
        if self.objective_direction == 'maximize':

            reaction_flux = simulate_reaction(self.__id, self.model, reaction, True)

            self.root.reactions_list = [(reaction, reaction.id, reaction_flux,
                                         list_reactants_ids(self.model, reaction.id),
                                         list_products_ids(self.model, reaction.id),
                                         get_reactants(self.model, reaction.id),
                                         get_products(self.model, reaction.id))]

        else:

            reaction_flux = simulate_reaction(self.__id, self.model, reaction, False)

            self.root.reactions_list = [(reaction, reaction.id, reaction_flux,
                                         list_products_ids(self.model, reaction.id),
                                         list_reactants_ids(self.model, reaction.id),
                                         get_products(self.model, reaction.id),
                                         get_reactants(self.model, reaction.id))]

        self.root.analysis = reaction_flux

    def set_reactions_list_to_node(self, node, last_reaction_list, reactant):

        """Creates the node reactions_list
        It populates the property reactions_list of the object node which is given as input
        with the reactions associated to the metabolite in the model
        the previous_reactions_list argument controls the appearances of repetitions.
        For instance, the metabolite protein should have only one reactions_list, the R_Protein.
        To do so, the R_Biomass should not be considered, even if the protein metabolite
        takes part into this reaction."""

        # the reactions_list of this metabolite should not contain the previous reactions_list Besides,
        # if the metabolite is a reactant, only reactions where the metabolite takes part as product should be
        # selected for reactions_list Otherwise, reactions where the metabolite takes part as reactant should be
        # selected for reactions_list

        if self.fast:

            reactions = get_reactions_by_role_fast(self.__id,
                                                   self.model,
                                                   node.id,
                                                   isReactant=reactant,
                                                   previous_reactions_list=last_reaction_list)

            node.reactions_list, node.other_reactions_list = reactions

        else:

            reactions = get_reactions_by_role(self.__id,
                                              self.model,
                                              node.id,
                                              isReactant=reactant,
                                              previous_reactions_list=last_reaction_list)

            node.reactions_list, node.other_reactions_list = reactions

    def create_next_nodes(self, node, leaf=False):

        if self.fast:

            if len(node.reactions_list) >= 20:

                node.isLeaf = True

            else:

                # let's iterate the property reactions_list of the node
                for reaction in node.reactions_list:
                    # for each reaction let's set the next nodes, namely the reactants of each reaction in
                    # reactions_list
                    self.create_next_nodes_by_reaction(node, reaction, leaf)

        else:

            # let's iterate the property reactions_list of the node
            for reaction in node.reactions_list:
                # for each reaction let's set the next nodes, namely the reactants of each reaction in reactions_list
                self.create_next_nodes_by_reaction(node, reaction, leaf)

    def create_next_nodes_by_reaction(self, node, reaction, leaf):

        reactants = reaction[5]
        products = reaction[6]

        for reactant in reactants:

            _next_node_hash = (reactant.id, reactant.name, reactant.compartment, 'reactant')
            # nextNode = node.has_next(reactant.id)
            nextNode = node.has_next_by_hash(_next_node_hash)

            if nextNode is None:
                # but if a metabolite with an equal name already exists,
                # it means that the metabolite is on a different compartment
                # let's create a composed name

                # nextNode = node.has_next_by_name(reactant.name)
                #
                # if nextNode is None:
                #     composed_name = reactant.name
                #
                # else:
                #     composed_name = reactant.name + ' ' + reactant.compartment

                self.build_node(identifier=reactant.id,
                                name=reactant.name,
                                compartment=reactant.compartment,
                                is_reactant=True,
                                previous_node=node,
                                leaf=leaf,
                                reactants=reactants,
                                products=products)

            # else:
            #
            #     if not nextNode.is_reactant:
            #         composed_name = reactant.name + ' ' + reactant.compartment + ' as reactant'
            #
            #         self.build_node(identifier=reactant.id,
            #                         name=composed_name,
            #                         compartment=reactant.compartment,
            #                         is_reactant=True,
            #                         previous_node=node,
            #                         leaf=leaf,
            #                         reactants=reactants,
            #                         products=products)

        # products version
        for product in products:

            _next_node_hash = (product.id, product.name, product.compartment, 'product')
            # isNode, nextNode = node.has_next(product.id)
            nextNode = node.has_next_by_hash(_next_node_hash)

            if nextNode is None:
                # but if a metabolite with an equal name already exists,
                # it means that the metabolite is on a different compartment
                # let's create a composed name

                # nextNode = node.has_next_by_name(product.name)
                #
                # if nextNode is None:
                #     composed_name = product.name
                #
                # else:
                #     composed_name = product.name + ' ' + product.compartment

                self.build_node(identifier=product.id,
                                name=product.name,
                                compartment=product.compartment,
                                is_reactant=False,
                                previous_node=node,
                                leaf=leaf,
                                reactants=reactants,
                                products=products)

            # else:
            #
            #     if nextNode.is_reactant:
            #         composed_name = product.name + ' ' + product.compartment + ' as product'
            #
            #         self.build_node(identifier=product.id,
            #                         name=composed_name,
            #                         compartment=product.compartment,
            #                         is_reactant=False,
            #                         previous_node=node,
            #                         leaf=leaf,
            #                         reactants=reactants,
            #                         products=products)

    def build_node(self, identifier, name, compartment, is_reactant, previous_node, leaf, reactants, products):

        # for each reactant in those reactions
        # create a new node using the id and name

        # the previous engine is similar to the next
        # the previous node of the next nodes is the current node itself

        # flag to control leafs
        # by default nodes are not leafs (see create_next_nodes)
        # populate tree controls the leafs

        # append next nodes of the current node, namely the reactants of the reactions where the node takes
        # part into
        # the products in the reactions_list reaction are the ones that might impair the flux

        next_node = Node(identifier=identifier,
                         name=name,
                         compartment=compartment,
                         is_reactant=is_reactant)
        next_node.previous = [previous_node]
        next_node.isLeaf = leaf
        previous_reactions_list = previous_node.reactions_list
        self.set_reactions_list_to_node(next_node, last_reaction_list=previous_reactions_list, reactant=is_reactant)
        previous_node.next.append(next_node)

        if is_reactant:
            next_node.analysis = simulate_reactants(self.__id, self.model, next_node, reactants, products)
        else:
            next_node.analysis = simulate_products(self.__id, self.model, next_node, reactants, products)

        return next_node

    def run(self, levels, fast=False):

        self.fast = fast

        if not self.root:
            self.set_root()

        if self.results is not None:
            self.results = None

        return self.populate_tree(levels)

    @timeout
    def populate_tree(self, levels):

        """
        create nodes associated with the root (biomass metabolite), namely each macromolecule or subunit,
        and subsequently
        """

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

            self.__populate_tree([self.root], level=0)

    def __populate_tree(self, nodes, level):

        """Recursive hidden method to populate the tree after the level 1
        After level 1 this function is called with the next nodes of root, namely the first next nodes
        for each node in the list of next nodes,
        the recursive function goes deep until it reaches the maximum level (controlled by level == self.levels)
        going through each level in a given node,
        the recursive function creates the next_nodes of this node and then get these for next recursion
        in the end, the final nodes are set as leafs
        """

        if level == self.levels:

            # exit for recursive function populate_sub_tree

            # each next node in the last level is set as leaf

            for node in nodes:
                node.isLeaf = True

            # finally, the recursion is ended
            return

        elif len(nodes) == 0:

            # exit for recursive function populate_sub_tree

            # each node without next node is set as leaf

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

                self.__populate_tree(nodes, level + 1)

    def get_tree(self):

        if self.results:
            return self.results

        # create empty tree to fill
        tree = {}

        # fill in tree starting with roots (those with no parent)
        self.__get_tree(tree, [self.root])

        self.results = tree

        return self.results

    def __get_tree(self, tree, nodes):

        # build a subtree for each child
        for child in nodes:

            _id = child.id
            _name = child.name
            _compartment = child.compartment
            _analysis = child.analysis
            _role = evaluate_side(child.is_reactant)
            _reactions_list = [(rxn[1], rxn[2], rxn[3], rxn[4])
                               for rxn in child.reactions_list]
            _other_reactions_list = [(rxn[1], rxn[2], rxn[3], rxn[4])
                                     for rxn in child.other_reactions_list]

            if _analysis:

                if child.has_next_nodes():

                    false_children = [next_child for next_child in child.get_next()
                                      if not next_child.analysis]

                    if len(false_children) > 0:
                        _analysis = False

            if not _analysis:

                true_rxns = [rxn[0] for rxn in _reactions_list if rxn[1]]

                if len(true_rxns) > 0:
                    _analysis = True

                if child.has_next_nodes() and not _analysis:

                    true_children = [next_child for next_child in child.get_next()
                                     if next_child.analysis]

                    if len(true_children) == len(child.get_next()):
                        _analysis = True

            # start new subtree
            tree[child.get_hash(stringify=True)] = {'identifier': _id,
                                                    'name': _name,
                                                    'compartment': _compartment,
                                                    'analysis': _analysis,
                                                    'role': _role,
                                                    'reactions': _reactions_list,
                                                    'other_reactions': _other_reactions_list,
                                                    'next': {}}

            # call recursively to build a subtree for current node
            self.__get_tree(tree[child.get_hash(stringify=True)]['next'], child.get_next())

    def write_results(self, results_f_name):

        results = self.get_tree()

        with open(results_f_name, "w") as jsonfile:
            json.dump(results, jsonfile)
