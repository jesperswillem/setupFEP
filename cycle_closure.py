# Name: Cycle Closure Error Estimator
"""
$Revision: 0.0 $

Compute the cycle closure errors, the minimum error pathway between all nodes
and nodes with experiment dG value (reference), the free energy associated with the
minimum error pathway for each ligand.

Copyright Schrodinger, LLC. All rights reserved.
"""
from __future__ import print_function

__all__ = [
    'CycleClosure',
]

import re
import csv
import hashlib
import Queue
import numpy as np
from math import sqrt
from itertools import chain
import networkx as nx

from schrodinger.application.scisol.packages.leadoptmap import tree

from . import graph
from .graph import Graph

# TODO: use edge properties
_EXP_DG_KEY = 'exp_dG'
_EXP_DDG_KEY = "exp_ddG"
_BENNETT_DDG_KEY = "bennett_ddG"
_BENNETT_ERROR_KEY = "bennett_error"
_CCC_DDG_KEY = "ccc_ddG"
_CCC_ERROR_KEY = "ccc_error"
_REL_DG_PATH_KEY = 'rel_dG_path'
_PRED_DG_KEY = 'pred_dG'
_PRED_DG_ERROR_KEY = 'pred_dG_error'
_GOOD_CYCLES_KEY = 'good_cycles'
_BAD_EDGE_KEY = "bad_edge"
_BAD_CYCLES_KEY = 'bad_cycles'
_LABEL_COUNT = 7

# GLOBALS
# This is the number of nodes we will span when searching for paths
# The larger this number the more cycles will be generated and the
# longer cycle closure will take to run
_NODE_CUTOFF = 6
_ACCEPTABLE_ERROR = 0.8
_EXPERIMENT_DG_ERROR_SQARE = 0.16  # 0.4 ^ 2

# TODO: do we need to expose this?
# This function was taken from leadoptmap.graph


def analyze(g, **kwargs):
    """
    Run cycle closure and update graph `g`

    :type g: `Graph`
    """
    cc = CycleClosure(g, **kwargs)
    cc.run_analysis()


class SubgraphCycleClosureAnalysis(object):
    """
    This is the core algorithm
    Uses a simplified graph with just weights and errors along edges
    Graph object is not modified
    """

    def __init__(self, graph, weighted_edges, bennett_error):
        """
        """
        self.graph = graph
        self.bennett_error = bennett_error
        self._hysteresis = []
        self._edges = weighted_edges
        self._cccm_error = [-1.0] * len(self._edges)
        self._best_paths = {}

    def process(self, rmse=True):
        """
        """
        path_list, base_path_list = self.find_close_cycles()
        cccm = self.populate_cccm(path_list)
        base_cccm = self.populate_cccm(base_path_list)
        ind_subgraphs = self.decomp_into_ind_subgraphs(base_cccm)
        independent_cccm_subgraph = self.getC3MForIndSubgraphs(
            ind_subgraphs, base_cccm)
        self.free_energy = self.solveCycleClosure(independent_cccm_subgraph,
                                                  ind_subgraphs)
        self.computeErrorsForCycles(cccm)
        if rmse:
            self.findBestPath()

    def path_dg_rmse(self, a_path, free_energy_dict):
        """
        """
        dg, path_sum = 0.0, 0.0
        for idx in range(len(a_path) - 1):
            A, B = a_path[idx], a_path[idx + 1]
            dg += free_energy_dict[(A, B)]
            path_sum += self.graph[A][B]['cccm_error2']
        return dg, path_sum

    def findBestPath(self):
        """
        First calculate estimated dG error for reference nodes (i.e., nodes
        with EXP_DG value); then find the best paths from reference nodes
        to all nodes in the graph
        """
        free_energy_dict = {}
        for idx, (n0, n1, _) in enumerate(self._edges):
            free_energy_dict[(n0, n1)] = self.free_energy[idx]
            free_energy_dict[(n1, n0)] = -1 * self.free_energy[idx]

        n = len(self.graph)  # number_of_nodes
        reference_nodes = [
            node for node, data in self.graph.nodes_iter(data=True)
            if (data.get(_EXP_DG_KEY) is not None and not node.is_ccc_excluded)
        ]
        pq = Queue.PriorityQueue()
        k = len(reference_nodes)
        if k != 1:
            paths = nx.shortest_path(self.graph, weight='cccm_error2')
        for i in range(k):
            node_i = reference_nodes[i]
            dg, sq_error = 0.0, 0.0
            for j in range(k):
                node = reference_nodes[j]
                dg += self.graph.node[node][_EXP_DG_KEY]
                sq_error += _EXPERIMENT_DG_ERROR_SQARE
                if i != j:
                    ddg, err = self.path_dg_rmse(paths[node_i][node],
                                                 free_energy_dict)
                    dg -= ddg
                    sq_error += err
            pq.put((sq_error / k**2, dg / k, [node_i]))

        # multiple target shortest path, using RMSE as weight
        # weighted BFS
        while n:
            error, dG_sum, path = pq.get()
            node = path[-1]
            if node not in self._best_paths:
                self._best_paths[node] = {
                    "error": sqrt(error),
                    "dG": dG_sum,
                    "path": path
                }
                n -= 1
                for neighbor in nx.all_neighbors(self.graph, node):
                    if neighbor not in self._best_paths:
                        edge = self.graph.get_edge(node, neighbor)
                        pq.put((error + edge.get_data('cccm_error2'),
                                dG_sum + free_energy_dict[(node, neighbor)],
                                path + [neighbor]))

    def computeErrorsForCycles(self, _cccm):
        """
        """

        # Each row is the length of edges so we can use the index where it is
        # non-zero
        # Add elements to the hysteresis list if they are unique.
        def rev(dictionary):
            '''
            Reverse the dictionary
            '''
            return {v: k for k, v in dictionary.iteritems()}

        hysteresis = {}
        for row in _cccm:
            lig_cycle = {}
            for idx, val in enumerate(row):
                if val != 0:
                    n1 = self._edges[idx][0]
                    n2 = self._edges[idx][1]
                    if val > 0:
                        lig_cycle[n1] = n2
                    else:
                        lig_cycle[n2] = n1
            cyc = frozenset(lig_cycle.iteritems())
            rev_cyc = frozenset(rev(lig_cycle).iteritems())
            if cyc not in hysteresis and rev_cyc not in hysteresis:
                e_list = [
                    val * edge[2] for val, edge in zip(row, self._edges) if val
                ]
                e_sum = abs(sum(e_list))
                e = e_sum / sqrt(len(e_list))

                hysteresis[cyc] = e_sum
                for idx, val in enumerate(row):
                    if val != 0 and e > self._cccm_error[idx]:
                        self._cccm_error[idx] = e

        for idx, (n0, n1, _) in enumerate(self._edges):
            e = self._cccm_error[idx] = max(self._cccm_error[idx],
                                            self.bennett_error[(n0, n1)])
            self.graph[n0][n1]['cccm_error2'] = e * e
        """
        Convert the cycle in self._hysteresis from dict to list
        sort _hysteresis by error
        """
        for cyc_set, err in hysteresis.iteritems():
            cyc_dict = dict(cyc_set)
            lig_cycle_lst = cyc_dict.keys()
            for i in range(len(cyc_set) - 1):
                lig_cycle_lst[i + 1] = cyc_dict[lig_cycle_lst[i]]
            self._hysteresis.append((lig_cycle_lst, err))
        self._hysteresis.sort(
            key=lambda x: x[1] / sqrt(len(x[0])), reverse=True)

    def find_close_cycles(self):
        """
        This will quickly find all simple cycles up to length
        _NODE_CUTOFF and base cycles in a graph.
        NetworkX is used to get the result
        """
        path_list = [
            path[:-1]
            for node in self.graph.nodes() for path in nx.all_simple_paths(
                self.graph, node, node, _NODE_CUTOFF) if len(path) > 3
        ]
        base_path_list = nx.cycle_basis(self.graph)
        return path_list, base_path_list

    def populate_cccm(self, path_list):
        """
        Will assign 1 and -1 for the flow in a cycle.
        the 11 edge example results in the following when the function returns
                           [[0,  1, -1, 0,  0, -1, -1, 0,  0,  0,  0],
                            [0,  0,  0, 0,  1, -1, -1, 0,  0,  0,  0],
                            [0,  0,  0, 0, -1,  1,  1, 0,  0,  0,  0],
                            [0,  0,  0, 0,  0,  0,  0, 0,  1,  1, -1],
                            [0,  0,  0, 0,  0,  0,  0, 0, -1, -1,  1],
                            [1,  0, -1, 1,  0, -1, -1, 0,  0,  0,  0],
                            [0,  1, -1, 0, -1,  0,  0, 0,  0,  0,  0],
                            [1,  0, -1, 1, -1,  0,  0, 0,  0,  0,  0],
                            [1, -1,  0, 1,  0,  0,  0, 0,  0,  0,  0]]
        """
        self.just_edges = {(e0, e1): i
                           for i, (e0, e1, _) in enumerate(self._edges)}

        def path2row(path):
            cccm_row = [0] * len(self._edges)
            l = len(path)
            for i, node in enumerate(path):
                edge = (node, path[(i + 1) % l])
                try:
                    idx = self.just_edges[edge]
                    cccm_row[idx] = 1
                except KeyError:
                    idx = self.just_edges[edge[::-1]]
                    cccm_row[idx] = -1
            return cccm_row

        return map(path2row, path_list)

    def decomp_into_ind_subgraphs(self, _base_cccm):
        """
        In this step we decompose our large graph into subgraphs. The algorithm
        that is used for this decomposition is we first go over all rows in the
        cccm and make a set for each row. In each set we record the column
        number if it is a non-zero element in the cccm. So for a 4x4 matrix we
        should end up with 4 sets. We then squash these sets. If a set has an
        intersection with some subgraphs, it is merged with all those subgraphs;
        otherwise it is added as it's own subgraph. At the end we
        convert the set to a list. Also if an edge element is not in any of the
        lists we add it to the result as itself.

        :param cccm: The cycle closure matrix represented by the following data
                     structure for the 11 edge example:
                     [[0,  1, -1, 0,  0, -1, -1, 0,  0,  0,  0],
                      [0,  0,  0, 0,  1, -1, -1, 0,  0,  0,  0],
                      [0,  0,  0, 0, -1,  1,  1, 0,  0,  0,  0],
                      [0,  0,  0, 0,  0,  0,  0, 0,  1,  1, -1],
                      [0,  0,  0, 0,  0,  0,  0, 0, -1, -1,  1],
                      [1,  0, -1, 1,  0, -1, -1, 0,  0,  0,  0],
                      [0,  1, -1, 0, -1,  0,  0, 0,  0,  0,  0],
                      [1,  0, -1, 1, -1,  0,  0, 0,  0,  0,  0],
                      [1, -1,  0, 1,  0,  0,  0, 0,  0,  0,  0]]
        :type cccm: A list of lists
        :param edges: A list of edges
        :type edges: list
        :return : A list of matrices each which represents an independent subgraph
        :rtype : A list of lists of lists
        """
        list_to_squash = [
            set(idx for idx, elem in enumerate(row) if elem != 0)
            for row in _base_cccm
        ]

        temp_list = []
        for current_set in list_to_squash:
            # isdisjoint should be faster than intersection
            merge_set = [
                i for i in range(len(temp_list))
                if not current_set.isdisjoint(temp_list[i])
            ]
            if len(merge_set) == 0:
                temp_list.append(current_set)
            else:
                j = merge_set[0]
                temp_list[j].update(current_set)
                for i in range(len(merge_set) - 1, 0, -1):
                    temp_list[j].update(temp_list.pop(merge_set[i]))
        temp_list = map(list, temp_list)

        col_set = set(range(len(self._edges)))
        diff_set = col_set.difference(*temp_list)
        temp_list.extend([[x] for x in diff_set])
        return temp_list

    def getC3MForIndSubgraphs(self, temp_list, _base_cccm):
        """
        In this step we get the C3M for each indpendent subgraph Ms.
        :param cccm: The cycle closure matrix represented by the following data
                     structure for the 11 edge example:
                     [[0,  1, -1, 0,  0, -1, -1, 0,  0,  0,  0],
                      [0,  0,  0, 0,  1, -1, -1, 0,  0,  0,  0],
                      [0,  0,  0, 0, -1,  1,  1, 0,  0,  0,  0],
                      [0,  0,  0, 0,  0,  0,  0, 0,  1,  1, -1],
                      [0,  0,  0, 0,  0,  0,  0, 0, -1, -1,  1],
                      [1,  0, -1, 1,  0, -1, -1, 0,  0,  0,  0],
                      [0,  1, -1, 0, -1,  0,  0, 0,  0,  0,  0],
                      [1,  0, -1, 1, -1,  0,  0, 0,  0,  0,  0],
                      [1, -1,  0, 1,  0,  0,  0, 0,  0,  0,  0]]
        :type cccm: a list of lists
        :param temp_list: the list of independent subgraphs for the 11 edge example
                          the list is equal to the following when passed in :
                          [[0, 1, 2, 3, 4, 5, 6], [8, 9, 10], [7]]
        :type temp_list: list of lists
        :return: a list of reduced matrices
        :rtype: a list of lists of lists
        """
        cccm_subgraph = []

        for a_list in temp_list:
            subgraph = []
            for row in _base_cccm:
                a_temp_list = [row[index] for index in a_list]
                # If all elements are zero don't append it...
                if any(elem != 0 for elem in a_temp_list):
                    subgraph.append(a_temp_list)
            cccm_subgraph.append(subgraph)

        return cccm_subgraph

    def solveCycleClosure(self, independent_cccm_subgraph, ind_subgraphs):
        """
        :param independent_cccm_subgraph : The reduced cccm subgraph for the full
                                           graph.  For the 11 edge example we have:
                                           [[[1, -1,  0, 1,  0,  0,  0],
                                             [1,  0, -1, 1, -1,  0,  0],
                                             [1,  0, -1, 1,  0, -1, -1]],
                                            [[-1, -1, 1]],
                                            [[]]
                                            ]
        :type independent_cccm_subgraph : A list of lists of lists
        :param edges : The edges for the full graph
        :type edges : a list of tuples
        :param ind_subgraphs : a list of lists of lists of integers for our 11 edge
                               example this is equal to
                               [[0, 1, 2, 3, 4, 5, 6], [8, 9, 10], [7]]
        :type ind_subgraphs :  a list of lists of lists
        :return : a list of free energies
        :rtype : a list of floats
        """
        # Init the F_array with the original computed delta delta G
        F_array = np.array([e[2] for e in self._edges], dtype=float)
        # For each list in ind_subgraphs
        # We don't want to include single edge graphs that cause Singular
        # Matrices
        ind_subgraphs_to_process = filter(lambda a_list: len(a_list) > 1,
                                          ind_subgraphs)
        for idx, a_list in enumerate(ind_subgraphs_to_process):
            # numpy rows
            num_subgraphs = len(independent_cccm_subgraph[idx])
            # numpy columns
            len_subgraph = len(independent_cccm_subgraph[idx][0])
            # Init our numpy matrix
            a_np_matrix = np.zeros(
                shape=(num_subgraphs, len_subgraph), dtype=float)
            neg_a_np_matrix = np.zeros(
                shape=(num_subgraphs, len_subgraph), dtype=float)
            E_array = np.array(range(len_subgraph), float)
            for list_idx in range(num_subgraphs):
                a_np_matrix[list_idx] = independent_cccm_subgraph[idx][list_idx]
                neg_a_np_matrix[list_idx] = \
                    [-1 * x for x in independent_cccm_subgraph[idx][list_idx]]
            for list_idx, e_idx in enumerate(a_list):
                E_array[list_idx] = self._edges[e_idx][2]

            dot_M_MT = np.dot(a_np_matrix, np.transpose(a_np_matrix))
            C_array = np.dot(neg_a_np_matrix, E_array)
            try:
                C_array = np.linalg.solve(dot_M_MT, C_array)
            except np.linalg.linalg.LinAlgError:
                msg = ("Singular matrix detected, cannot complete free energy"
                       " calculation. Please check the ddG input values.")
                print(msg)
                return F_array

            result = E_array + np.dot(np.transpose(a_np_matrix), C_array)
            for f_idx, free_energy in enumerate(result):
                F_array[a_list[f_idx]] = free_energy
        return F_array


class CycleClosure(object):
    """
    Interface to cycle closure.
    """

    def __init__(self, graph=None, results_fn=None, results_fmp=None):
        """
        Initialize a new CycleClosure object.

        :param graph: graph or None
        :type graph: `Graph` or None
        :param results_fn: Filename for output or None for no output
        :type results_fn: string or None
        """
        self.graph = graph
        self.results_fn = results_fn
        self.results_fmp = results_fmp

    def readCsv(self, filename):
        """
        This will read the input from a csv file which is prepared by the user.
        We build a graph which uses the hash of the titles as the node id's just
        as a normal run through the mapper would generate.  We also populate the
        edge with the same key that mapper would so we can use the same framework
        to compute the ccc_ddg and ccc_error.  We must take care when generating
        the key and make sure the sign obbeys hexA < hexB.

        :param filename: The name of the csv file to read the data from
        :type filename: string
        """

        mol_dict_map = {}
        rev_mol_dict_map = {}
        csv_reader = csv.reader(open(filename))
        G = Graph()
        for idx, row in enumerate(csv_reader):
            if idx == 0:
                check_csv_header(row)
            else:
                try:
                    mol1, mol2, ddG = row[0].strip().replace(" ", "_"), \
                        row[1].strip().replace(" ", "_"), \
                        row[2].strip()
                except IndexError:
                    pass
                else:
                    hexA = hashlib.sha1(mol1).hexdigest()
                    hexB = hashlib.sha1(mol2).hexdigest()
                    sign = 1
                    if hexA > hexB:
                        # We enforce the standard that the perturbation
                        # should be hexA -> hexB by hexA < hexB
                        sign = -1

                    if mol1 not in mol_dict_map:
                        node = graph.Node(G, id=hexA)
                        mol_dict_map[mol1] = node
                        rev_mol_dict_map[node] = mol1
                        G.add_node(node)

                    if mol2 not in mol_dict_map:
                        node = graph.Node(G, id=hexB)
                        mol_dict_map[mol2] = node
                        rev_mol_dict_map[node] = mol2
                        G.add_node(node)

                    try:
                        e = row[3].strip()
                    except IndexError:
                        e = 0.0

                    G.add_edge(
                        mol_dict_map[mol1],
                        mol_dict_map[mol2],
                        bennett_ddG=sign * float(ddG),
                        bennett_error=float(e))

        for (node, prop) in G.nodes(data=True):
            prop['title'] = rev_mol_dict_map[node]
            prop['label'] = short_label(node.id)
        for (n0, n1), data in G.edges_iter(data=True):
            label = short_label(n0.id) + "-" + short_label(n1.id)
            data['label'] = label
        self.graph = G

    def run_analysis(self):
        """
        Run Cycle Closure and output results to results_fn if set.
        """
        self.analysis = CycleClosureAnalysis(self.graph)
        self.analysis.populateGraph()
        self.analysis.processSubgraphs()
        self.analysis.collectResults()
        if self.results_fn:
            self.analysis.writeResults(self.results_fn)
        if self.results_fmp:
            self.analysis.graph.write(self.results_fmp)


class CycleClosureAnalysis(object):
    """
    """

    def __init__(self, graph):
        """
        """
        self.graph = graph
        self.map, self.rev_map = get_simple_map(self.graph)
        self.weighted_edges = []
        self.bennett_error = {}
        # This will contain lists of nodes in each disconnected graph
        self._processed_subgraphs = []

    def populateGraph(self):
        """
        Populate weighted edge and bennett error
        """
        _, self.rev_map = get_simple_map(self.graph)
        for edge in self.graph.edges_iter():
            n0, n1 = edge
            try:
                ddg = edge.bennett_ddg.val
                err = edge.bennett_ddg.unc
            except AttributeError:
                t1 = self._get_node_attribute(self.graph, n0.id, 'title')
                t2 = self._get_node_attribute(self.graph, n1.id, 'title')
                print("WARNING: Edge %s-%s, does not have a ddG value!" % (t1,
                                                                           t2))
            else:
                self.weighted_edges.append((n0, n1, ddg))
                self.bennett_error[(n0, n1)] = self.bennett_error[(n1,
                                                                   n0)] = err

    def processSubgraphs(self):
        """
        Run cycle closure on disconnected subgraphs
        """
        # generate a graph that all edges has ddG value
        temp_graph = self._create_graph_from_weighted_edges(self.weighted_edges)

        for nodes in nx.connected_components(temp_graph):
            weighted_edges = [
                w for w in self.weighted_edges
                if w[0] in nodes and w[1] in nodes
            ]
            subgraph = self._create_graph_from_weighted_edges(weighted_edges)

            no_reference = True
            for n, d in subgraph.nodes_iter(data=True):
                exp_dg = n.exp_dg
                if exp_dg is not None:
                    d[_EXP_DG_KEY] = exp_dg
                    no_reference = False
            if no_reference:
                root = get_root_node(subgraph)
                subgraph.node[root][_EXP_DG_KEY] = 0.0

            cc = SubgraphCycleClosureAnalysis(subgraph, weighted_edges,
                                              self.bennett_error)
            cc.process()
            self._processed_subgraphs.append(cc)

    def collectResults(self):
        """
        First clear previous properties in graph, then
        add computed properties to graph
        """

        def try_del(d, k):
            d.pop(k, None)

        for _, d in self.graph.edges_iter(data=True):
            try_del(d, _CCC_ERROR_KEY)
            try_del(d, _CCC_DDG_KEY)
            try_del(d, _BENNETT_ERROR_KEY)
            try_del(d, _BENNETT_DDG_KEY)
            try_del(d, _EXP_DDG_KEY)
            try_del(d, _BAD_EDGE_KEY)

        for _, d in self.graph.nodes_iter(data=True):
            try_del(d, _PRED_DG_ERROR_KEY)
            try_del(d, _REL_DG_PATH_KEY)
            try_del(d, _PRED_DG_KEY)

        self.graph.graph[_GOOD_CYCLES_KEY] = {}
        self.graph.graph[_BAD_CYCLES_KEY] = {}

        for cc in self._processed_subgraphs:
            # Adds to full graph not the subgraph
            self.addEdgePropertiesToGraph(cc)
            self.addNodePropertiesToGraph(cc)
            self.addGraphPropertiesToGraph(cc)
            self.label_bad_edges(cc)

    def _get_node_attribute(self, graph, node_id, attribute):
        """
        """
        node = self.rev_map[node_id]
        try:
            return graph.node[node][attribute]
        except KeyError:
            if attribute == 'title':
                return graph.node[node]['mmct'].title

    def writeResults(self, filename):
        """
        """
        with open(filename, 'w') as fh:
            for cc in self._processed_subgraphs:
                fh.write(
                    "\nLigand_1, Ligand_2, ddG(kcal/mol), Cycle Closure Error "
                    "Along Edge, Input ddG(kcal/mol), Input Bennett error, "
                    "Hexcode Lig1, Hexcode Lig2 \n")
                fh.write(
                    "---------------------------------------------------\n")
                for l in self.get_formatted_ddG_results(cc):
                    fh.write(
                        "%5s, %5s, %4.2f, %4.2f, %4.2f, %4.2f, %s, %s\n" % l)

                fh.write("\nHysteresis along cycle, Cycle\n")
                fh.write(
                    "---------------------------------------------------\n")
                for cyc, hys_error in cc._hysteresis:

                    a_new_list = [
                        self._get_node_attribute(self.graph, n.id, 'title')
                        for n in cyc
                    ]
                    if hys_error / sqrt(len(cyc)) > _ACCEPTABLE_ERROR:
                        fh.write("%4.2f, [%s] *\n" % (hys_error,
                                                      ' '.join(a_new_list)))
                    else:
                        fh.write("%4.2f, [%s]\n" % (hys_error,
                                                    ' '.join(a_new_list)))

                fh.write(
                    "* Large hysteresis error along edge has been detected, "
                    "please check your perturbations within the cycle.\n")

                fh.write("\nLigand, ddG(kcal/mol), Error along path, Path\n")
                fh.write(
                    "---------------------------------------------------\n")
                for l in self.get_formatted_dG_results(cc):
                    fh.write("%5s, %4.2f, %4.2f, %s\n" % l)

    def get_formatted_dG_results(self, cc):
        """
        Returns a formatted sorted list of dG results for each node, with the relative path

        :return: List of dG values for each node
                [
                [ node, dG, dG_err, path ],
                ...
            ]
        :rtype: list of list of dG results per node
        """
        to_sort = []
        for node in cc.graph.nodes():
            n = self.graph.node[node]
            title = self._get_node_attribute(self.graph,
                                             n[_REL_DG_PATH_KEY][-1], 'title')
            path = '--->'.join([
                self._get_node_attribute(self.graph, a, 'title')
                for a in n[_REL_DG_PATH_KEY]
            ])
            to_sort.append((title, n[_PRED_DG_KEY], n[_PRED_DG_ERROR_KEY],
                            path))
        to_sort.sort(key=lambda tup: tup[1])
        return to_sort

    def get_formatted_ddG_results(self, cc):
        """
        Returns a formatted list of ddG results for each edge

        :return: List of ddG values for each edge
                [
                    [ node_from, node_to, ccc_ddG, ccc_err, bennett_ddG, bennett_error,
                      hexid0, hexid1 ],
                    ...
                ]
        :rtype: list of list of ddG results per edge
        """
        results = []
        for edge in cc.graph.edges():
            n0, n1 = edge
            edge_dict = self.graph.get_edge_data(n0, n1)
            edge_property = (self._get_node_attribute(
                self.graph, n0.id, 'title'), self._get_node_attribute(
                    self.graph, n1.id, 'title'), edge_dict[_CCC_DDG_KEY],
                             edge_dict[_CCC_ERROR_KEY],
                             edge_dict[_BENNETT_DDG_KEY],
                             edge_dict[_BENNETT_ERROR_KEY], short_label(n0.id),
                             short_label(n1.id))
            results.append(edge_property)
        return results

    def _create_graph_from_weighted_edges(self, weighted_edges):
        g = Graph()
        all_nodes = set(
            list(chain.from_iterable((u, v) for (u, v, _) in weighted_edges)))
        all_node_pairs = set((u, v) for (u, v, _) in weighted_edges)

        for node in all_nodes:
            g.add_node(node)
        for pair in all_node_pairs:
            g.add_edge(*pair)
        g.add_weighted_edges_from(weighted_edges)
        return g

    def label_bad_edges(self, cc, count=1):
        """
        for the worst cycle (biggest hystersis), label edges with biggest difference between
        bennet_ddG and CycleClosure_ddG as potential bad cycle.
        Hypothetically delete these edges, and repeat.
        GUI will display this information and prompt user to consider delete them.

        :type cc: SubgraphCycleClosureAnalysis object
        :type count: int
        :return: None   (bad edges is labeled on the graph)
        """
        while cc._hysteresis:
            cyc, error = cc._hysteresis[0]
            length = len(cyc)
            if error / sqrt(length) > _ACCEPTABLE_ERROR:
                bad_edges = []
                for i in range(length):
                    n0 = cyc[i]
                    n1 = cyc[(i + 1) % length]
                    try:
                        idx = cc.just_edges[(n0, n1)]
                    except KeyError:
                        idx = cc.just_edges[(n1, n0)]
                    n0, n1, bennet_ddG = cc._edges[idx]
                    ccc_ddG = cc.free_energy[idx]
                    bad_edges.append((idx, n0, n1, abs(bennet_ddG - ccc_ddG)))
                # find the biggest ddG deviation value
                m = max(bad_edges, key=lambda x: x[3])[3]
                deleted_edges = []
                for idx, n0, n1, err in bad_edges:
                    # these are floating numbers, and test equality is not
                    # reliable
                    if err > m - 0.00000001:
                        self.graph[n0][n1][_BAD_EDGE_KEY] = count
                        deleted_edges.append(idx)
                _weighted_edges = [
                    e for i, e in enumerate(cc._edges) if i not in deleted_edges
                ]
                subgraph = self._create_graph_from_weighted_edges(
                    _weighted_edges)
                cc = SubgraphCycleClosureAnalysis(subgraph, _weighted_edges,
                                                  self.bennett_error)
                cc.process(False)
                count += 1
            else:
                break

    def addGraphPropertiesToGraph(self, cc):
        """
        :param cc : object which contains information to be added to graph
        :type cc :  SubgraphCycleClosureAnalysis object
        """
        good_cycles = {}
        bad_cycles = {}
        for cyc, hys_val in cc._hysteresis:
            hex_list = [w.id for w in cyc]
            cycle_nodes = tuple(hex_list)
            if hys_val / sqrt(len(cyc)) < _ACCEPTABLE_ERROR:
                good_cycles[cycle_nodes] = "%4.2f" % hys_val
            else:
                bad_cycles[cycle_nodes] = "%4.2f" % hys_val

        self.graph.graph[_GOOD_CYCLES_KEY].update(good_cycles)
        self.graph.graph[_BAD_CYCLES_KEY].update(bad_cycles)

    def addNodePropertiesToGraph(self, cc):
        """
        :param cc : object which contains information to be added to graph
        :type cc :  SubgraphCycleClosureAnalysis object
        """
        for node, path in cc._best_paths.iteritems():
            hex_path = [a.id for a in path['path']]
            # Use the same Key that Tom is using for GUI
            pred_dg = float("%4.2f" % path['dG'])
            self.addPropToGraphNode(node, _PRED_DG_KEY, pred_dg)
            error = float("%4.2f" % path['error'])
            self.addPropToGraphNode(node, _PRED_DG_ERROR_KEY, error)
            self.addPropToGraphNode(node, _REL_DG_PATH_KEY, hex_path)

    def addEdgePropertiesToGraph(self, cc):
        """
        :param cc : object which contains information to be added to graph
        :type cc :  SubgraphCycleClosureAnalysis object
        """
        for idx, (n0, n1, dg) in enumerate(cc._edges):
            ccc_free = float("%4.2f" % cc.free_energy[idx])
            ccc_err = float("%4.2f" % cc._cccm_error[idx])
            ben_free = float("%4.2f" % dg)
            ben_err = float("%4.2f" % self.bennett_error[(n0, n1)])

            # desomnd-4911
            if ccc_err < 0:
                ccc_err = float('nan')

            self.addPropToGraphEdge(n0, n1, _CCC_ERROR_KEY, ccc_err)
            self.addPropToGraphEdge(n0, n1, _CCC_DDG_KEY, ccc_free)
            self.addPropToGraphEdge(n0, n1, _BENNETT_ERROR_KEY, ben_err)
            self.addPropToGraphEdge(n0, n1, _BENNETT_DDG_KEY, ben_free)

        for (n0, n1), data in self.graph.edges_iter(data=True):
            expdg1 = n0.exp_dg
            expdg2 = n1.exp_dg
            if expdg1 is not None and expdg2 is not None:
                data[_EXP_DDG_KEY] = expdg2 - expdg1

    def addPropToGraphEdge(self, node0, node1, key, val):
        """
        """
        edge = self.graph.get_edge(node0, node1)
        edge.set_data(key, val)

    def addPropToGraphNode(self, node, key, val):
        """
        """
        self.graph.node[node][key] = val


def short_label(long_hex):
    """
    Will return the shortened version used throughout multisim
    :param long_hex : The long hexidecimal string
    :type long_hex : string
    """
    return long_hex[:_LABEL_COUNT]


def get_simple_map(graph):
    """
    This just maps the nodes full hex id to an integer value
    The integer values for the nodes is used in the algorithm
    for computing cc
    """
    node_map, rev_map = {}, {}
    for idx, node in enumerate(graph.nodes()):
        node_map[node] = idx
        rev_map[node.id] = node
    return node_map, rev_map


def get_root_node(tree_graph, base=0):
    """
    If we fail to find a node return first node as the rooth
    :param tree_graph: The tree to generate the graph for finding the root.
    :type tree_graph: `graph.Graph`
        return: `graph.Node`. Will return first `Node` if can not build the tree
        rtype: `graph.Node`
    """
    pickle_graph = tree_graph.copy()
    atree = tree.graph2tree(pickle_graph)
    for node in str(atree).split("\n"):
        node_num = re.sub('[+-]', '', node)
        for node2 in tree_graph.nodes():
            if node2.id.startswith(node_num):
                return node2
    return tree_graph.nodes()[0]


def check_csv_header(row):
    """
    This will check to make sure we have a header...the third value
    """
    try:
        ddG = row[2].strip()
    except IndexError:
        print("Could not read header in csv file")
    else:
        try:
            float(ddG)
        except ValueError:
            pass
        else:
            print("You must include a header in your csv file, for example\n "
                  "lig1, lig2, ddG, error")
