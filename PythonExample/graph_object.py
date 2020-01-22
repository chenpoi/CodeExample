
"""
This scripts is to remove redundancy in a collections of sequences.
input:
    all sequences ids in the collections
    self blast output in format 6 (blast agains it self)
output:
    a tsv file, which group each sequence into a cluster
"""

import re
import argparse
import operator
import time
import pickle
import math
import networkx


class Graph:
    """
    it's the base class of Graph object, containing two class as component, Vertex and Edge
    The structure of graph is based on Adjacency Map Structure.
    it also contains some fundamental methods
    """
    class Vertex:
        """
        _element: the name of vertex
        _visited: the status of vertex, not used in this case
        """
        __slots__ = '_element', '_visited'

        def __init__(self, x):
            self._element = x
            self._visited = False

        def element(self):
            return self._element

        def is_visit(self):
            return self._visited

        def get_visited(self):
            self._visited = True

        def __hash__(self):
            return hash(id(self))

    class Edge:
        __slots__ = '_origin', '_destination'

        def __init__(self, u, v):
            """
            :param u: the origin vertex of edge
            :param v: the end vertex of edge
            """
            self._origin = u
            self._destination = v

        def endpoints(self):
            return self._origin, self._destination

        def opposite(self, v):
            return self._destination if v is self._origin else self._origin

        def __hash__(self):
            return hash((self._origin, self._destination))

    def __init__(self, directed=False):
        """
        _outgoing is a dictionary of dictionary. the value inside is a certain edge.
            the outer key is the origin vertex and the inside key is the destination of edge
        _incoming is just opposite to _outgoing.
        _vertex_dict will give you the vertex object when you have the name of vertex
        :param directed: if it's a directed graph
        """
        self._outgoing = {}
        self._incoming = {} if directed else self._outgoing
        self._vertex_dict = {}

    def check_reverse_edge(self, e):
        try:
            self._outgoing[e.endpoints()[1]][e.endpoints()[0]]
            return True
        except KeyError:
            return False

    def is_directed(self):
        return self._incoming is not self._outgoing

    def vertices(self):
        return self._outgoing.keys()

    def edge_count(self):
        total = sum(len(self._outgoing[v]) for v in self._outgoing)
        return total if self.is_directed() else total // 2

    def edges(self):
        result = set()
        for secondaty_map in self._outgoing.values():
            result.update(secondaty_map.values())
        return result

    def get_edge(self, u, v):
        u_vertex = self._vertex_dict[u]
        v_vertex = self._vertex_dict[v]
        return self._outgoing[u_vertex].get(v_vertex)

    def degree(self, v, outgoing=True):
        v_vertex = self._vertex_dict[v]
        adj = self._outgoing if outgoing else self._incoming
        return len(adj[v_vertex])

    def incident_edges(self, v, outgoing=True):
        if isinstance(v, self.Vertex):
            v_vertex = v
        else:
            v_vertex = self._vertex_dict[v]
        adj = self._outgoing if outgoing else self._incoming
        for edge in adj[v_vertex].values():
            yield edge

    def incident_vertex(self, v, outgoing=True):
        return self._outgoing[v].keys() if outgoing else self._incoming[v].keys()

    def insert_vertex(self, x=None):
        v = self.Vertex(x)
        self._vertex_dict[x] = v
        self._outgoing[v] = {}
        if self.is_directed():
            self._incoming[v] = {}
        return v

    def insert_edge(self, u_element, v_element, **kwargs):
        u_vertex = self._vertex_dict[u_element]
        v_vertex = self._vertex_dict[v_element]
        e = self.Edge(u_vertex, v_vertex)
        self._outgoing[u_vertex][v_vertex] = e
        self._incoming[v_vertex][u_vertex] = e


class BlastGraph(Graph):
    """
    it's the graph designed for our project to store the result of blast.
    it contains a BlastEdge class, which is customized version of base Edge class. the slots '_pident' and '_qcovs' are
    from blast stats and '_score' is defined as the multiply of '_pident' and '_qcovs'
    this class also has three methods for constructing the graph from a particular format of blast result.
    """
    class BlastEdge(Graph.Edge):
        __slots__ = '_origin', '_destination', '_pident', '_qcovs', '_evalue', '_score'

        def __init__(self, u, v, pident, qcovs, evalue=None):
            self._pident = pident
            self._qcovs = qcovs
            self._evalue = evalue
            if evalue:
                if float(evalue) == 0:
                    evalue = "1e-180"
                self._score = 1 - (math.log10(float(evalue)) + 180) / (1 + 180)
                #  1 / (1 + math.exp(0.1 * (math.log10(float(evalue)) + 30)))
            else:
                self._score = float(pident) * float(qcovs) / 10000

            super().__init__(u, v)

        def get_pident(self):
            return self._pident

        def get_qcovs(self):
            return self._qcovs

        def get_score(self):
            return self._score

    def insert_edge(self, u, v, **kwargs):
        cutoff = min(kwargs["thresholds"])
        del kwargs["thresholds"]
        u_vertex = self._vertex_dict[u]
        v_vertex = self._vertex_dict[v]
        e = self.BlastEdge(u_vertex, v_vertex, **kwargs)
        if e.get_score() < cutoff:
            return 0
        else:
            self._outgoing[u_vertex][v_vertex] = e
            self._incoming[v_vertex][u_vertex] = e
            return e

    def read_vertex(self, file):
        with open(file, "r") as vertex_file:
            for cnt, line in enumerate(vertex_file):
                self.insert_vertex(line.rstrip())

    def read_edge(self, file, qseqid, sseqid, pident, qcovs, evalue, thresholds):
        with open(file, "r") as edge_file:
            for cnt, line in enumerate(edge_file):
                words = re.split(r"\s+", line.rstrip())
                query_id = words[qseqid]
                subject_id = words[sseqid]

                if query_id != subject_id:
                    # if evalue is specified, then enter evalue mode
                    if evalue:
                        self.insert_edge(query_id, subject_id,
                                         pident=words[pident], qcovs=words[qcovs],
                                         evalue=words[evalue], thresholds=thresholds)
                    else:
                        self.insert_edge(query_id, subject_id,
                                         pident=words[pident], qcovs=words[qcovs],
                                         thresholds=thresholds)

        '''
        for u in self._outgoing.keys():
            for v in self._outgoing[u].keys():
                if self.self._outgoing[u][v].get_score() == 1:
                    try:
                        if self.self._outgoing[v][u].get_score() != 1:
                            del self.self._outgoing[v][u]
                    except KeyError:
                        continue
        '''


# Graph Based Sequence Clustering Algorithm
class GBSCA(BlastGraph):
    """
    it's the main body of our Graph Based Sequence Clustering Algorithm(GBSCA). the class Cluster contain method of
    cluster constructor, checking if a node can be added to that cluster and checking if two clusters can be merged.
    the '_cluster_nodes' slots contain a dictionary of node and their outdegree. the '_path_distance' is the shortest
    path between any two node, if reachable, in the cluster.
    GBSCA class itself also has a method dispatch to guide the analysis into certain method. for the four dispatched
    method, they all consist of four parts: firstly, initial the tested nodes and clusters and do some simple test;
    Secondly, get the in_edges and out_edges between them. Thirdly, check it they can be added/merged. Fourth, update
    the cluster and node information in GBSCA.
    """
    class Cluster:
        __slots__ = "_cluster_nodes",  "_path_distance", "_cluster_representative", "_verbose"

        def __init__(self, u, v, e, verbose):
            self._cluster_nodes = {u: 1, v: 0}  # node as key and outdegree as value
            self._path_distance = {}
            self._cluster_representative = []

            self._path_distance[u] = {v: e.get_score()}
            self._path_distance[v] = {u: 0}
            self._verbose = verbose

        def set_verbose(self, verbose):
            self._verbose = verbose

        def get_cluster_nodes(self):
            return self._cluster_nodes

        def get_outdegree(self, u):
            return self._cluster_nodes[u]

        def get_path_distance(self, u):
            return self._path_distance[u]

        def get_shortest_distance(self, u, v=None, incoming=False):
            if incoming:
                return self._path_distance[v][u]
            else:
                return self._path_distance[u][v]

        def initialize_node_insertion(self, u):
            self._cluster_nodes[u] = 0
            self._path_distance[u] = {}

        def set_shortest_distance(self, u, v, dist):
            self._path_distance[u][v] = dist

        def get_reachable_nodes(self, u, incoming=False):
            """
            :param u: certain node in the cluster
            :param incoming: it's for direction. if you want to regard u as the destination of the path,
                                make "incoming" True.
            :return: all the node that can connect to u in certain direction
            """
            if incoming:
                return [v for v in self.get_cluster_nodes().keys() if u in self._path_distance[v].keys()]
            else:
                return self._path_distance[u].keys()

        def calculate_one_to_cluster_distance(self, v, bridge_dict, edges_check, incoming=False):
            """
            :param v: any vertex in this cluster
            :param bridge_dict: the collection of edges connect the nodes in this cluster and outside. its keys are the
                nodes in this clusters (so there is a possibility that v may appear in bridge_dict too). its values are
                the bridge edges.
            :param edges_check: threshold to check if there is any need to do a score multiply
            :param incoming: get the direction information. if the bridge edges are from outside into the clusters,
                incoming should be true. otherwise, it's false
            :return: the shortest distance from a bridges edge to node v inside the cluster.
            """
            if v in bridge_dict.keys():
                dist = bridge_dict[v].get_score()
            else:
                dist = max([bridge_dict[t].get_score() * self.get_shortest_distance(v, t, incoming)
                            for t in bridge_dict.keys()
                            if self.get_shortest_distance(v, t, incoming) > edges_check] + [0])
            return dist

        def check_and_insert_node(self, u, cluster_threshold, in_edges, out_edges, double_edges,
                                  in_edges_max_score, out_edges_max_score):
            """
            temp_in and temp_out will store all the shortest distance between u and any node in the clusters,
            if reachable. if temp_in and temp_out pass the cluster_threshold test, they can also be used to update
            "_path_distance" when u is added. At the same time, outdegree of any node in the clusters is also updated.
            :param u: the node to be checked
            :param cluster_threshold: threshold for the minimum _path_distance in cluster
            :param in_edges: all the edges from u to the clusters
            :param out_edges: all the edges from the clusters to u
            :param in_edges_max_score: threshold to check if there is any need to do a score multiply
            :param out_edges_max_score: threshold to check if there is any need to do a score multiply
            :return: if u is added to the cluster or not
            """
            in_edges_check = cluster_threshold / in_edges_max_score[0]
            out_edges_check = cluster_threshold / out_edges_max_score[0]

            temp_in_nodes = {e.endpoints()[1]: e for e in in_edges}
            temp_out_nodes = {e.endpoints()[0]: e for e in out_edges}
            temp_in = {}
            temp_out = {}

            for v in self.get_cluster_nodes().keys():
                dist_in = self.calculate_one_to_cluster_distance(
                    v, temp_in_nodes, in_edges_check, incoming=True
                )
                dist_out = self.calculate_one_to_cluster_distance(
                    v, temp_out_nodes, out_edges_check
                )
                if 0 < dist_in < cluster_threshold and 0 < dist_out < cluster_threshold:
                    if self._verbose:
                        print("first cluster is: " +
                              "; ".join([j.element() for j in self.get_cluster_nodes().keys()]) +
                              "\nThe failed insertion is from: " + u.element() + " to " + v.element() +
                              "\nforward distance is " + str(dist_in) + "; backward distance is " + str(dist_out))
                    # fail
                    return 0
                else:
                    temp_in.update({v: dist_in})
                    temp_out.update({v: dist_out})

            # success
            self.initialize_node_insertion(u)

            self._cluster_nodes[u] += len([1 for i in in_edges
                                           if i not in double_edges])
            for i in out_edges:
                if i not in double_edges:
                    self._cluster_nodes[i.endpoints()[0]] += 1

            for v, dist in temp_in.items():
                self._path_distance[u][v] = dist
            for v, dist in temp_out.items():
                self._path_distance[v][u] = dist
            return 1

        def check_and_merge_cluster(self, clu_v, cluster_threshold, in_edges, out_edges, double_edges,
                                    in_edges_max_score, out_edges_max_score):
            """
            the idea is similar to check_and_insert_node method. but expand the one-to-many relation to a many-to-many
            relation.
            the border nodes are the ones who are in both the endpoints of bridges edges and clu_v
            :param clu_v: the other cluster
            :param cluster_threshold: threshold for the minimum _path_distance in cluster
            :param in_edges: all the edges from clu_u to clu_v
            :param out_edges: all the edges from clu_v to clu_u
            :param in_edges_max_score: threshold to check if there is any need to do a score multiply
            :param out_edges_max_score: threshold to check if there is any need to do a score multiply
            :return: if two clusters are merged or not
            """

            in_edges_check = cluster_threshold / in_edges_max_score[0]
            out_edges_check = cluster_threshold / out_edges_max_score[0]

            border_nodes_set = set([e.endpoints()[0] for e in in_edges] +
                                   [e.endpoints()[1] for e in out_edges])

            temp_border_dist = {}
            for i in border_nodes_set:
                temp_in_nodes = {e.endpoints()[1]: e for e in in_edges if e.endpoints()[0] == i}
                temp_out_nodes = {e.endpoints()[0]: e for e in out_edges if e.endpoints()[1] == i}
                temp_in = {}
                temp_out = {}

                for v in clu_v.get_cluster_nodes().keys():
                    dist_in = clu_v.calculate_one_to_cluster_distance(
                        v, temp_in_nodes, in_edges_check, incoming=True
                    )
                    dist_out = clu_v.calculate_one_to_cluster_distance(
                        v, temp_out_nodes, out_edges_check
                    )
                    if 0 < dist_in < cluster_threshold and 0 < dist_out < cluster_threshold:
                        # fail
                        if self._verbose:
                            print("The failed insertion is from: " + i.element() + " to " + v.element())
                        return 0
                    else:
                        temp_in.update({(i, v): dist_in})
                        temp_out.update({(v, i): dist_out})

                temp_border_dist.update({i: [temp_in, temp_out]})

            temp_distance = {"in": {}, "out": {}}
            for i in self.get_cluster_nodes().keys():
                if i in border_nodes_set:
                    temp_distance["in"].update({
                        pair_key: dist for pair_key, dist in temp_border_dist[i][0].items()
                    })
                    temp_distance["out"].update({
                        pair_key: dist for pair_key, dist in temp_border_dist[i][1].items()
                    })
                else:
                    for v in clu_v.get_cluster_nodes().keys():
                        # only if the edge with score more than edge_check will get multiplication,
                        # which is the most time-consuming step.
                        dist_in = max([self.get_shortest_distance(i, border_node) *
                                       temp_border_dist[border_node][0][(border_node, v)]
                                       for border_node in border_nodes_set
                                       if self.get_shortest_distance(i, border_node) > in_edges_check] + [0])
                        dist_out = max([self.get_shortest_distance(border_node, i) *
                                        temp_border_dist[border_node][1][(v, border_node)]
                                        for border_node in border_nodes_set
                                        if self.get_shortest_distance(border_node, i) > out_edges_check] + [0])
                        if 0 < dist_in < cluster_threshold and 0 < dist_out < cluster_threshold:
                            if self._verbose:
                                print("The failed insertion is from: " + i.element() + " to " + v.element())
                            # fail
                            return 0
                        else:
                            temp_distance["in"].update({(i, v): dist_in})
                            temp_distance["out"].update({(v, i): dist_out})

            # success
            for i in clu_v.get_cluster_nodes().keys():
                self._cluster_nodes[i] = clu_v.get_outdegree(i)
                self._path_distance[i] = clu_v.get_path_distance(i)

            for pair_key in temp_distance["in"]:
                self._path_distance[pair_key[0]][pair_key[1]] = temp_distance["in"][pair_key]

            for pair_key in temp_distance["out"]:
                self._path_distance[pair_key[0]][pair_key[1]] = temp_distance["out"][pair_key]

            for i in out_edges:
                if i not in double_edges:
                    self._cluster_nodes[i.endpoints()[0]] += 1
            for i in in_edges:
                if i not in double_edges:
                    self._cluster_nodes[i.endpoints()[0]] += 1
            return 1

        def show_cluster_component(self, out_file=None):
            """
            This method will show you the nodes, edges information in a cluster
            :param out_file: if you want to push the output to a file, specify the file name here
            :return: None
            """
            nodes = self.get_cluster_nodes().keys()
            edges = []
            for i in nodes:
                edges.append("\n".join([i.element() + "\t" + j.element() + "\t" + str(self.get_shortest_distance(i, j))
                                        for j in self.get_reachable_nodes(i)]))

            msg = "nodes are\n" + "\n".join(map(operator.methodcaller('element'), nodes)) + \
                  "\nedges are\n" + "\n".join(edges)

            if out_file:
                out_fh = open(out_file, "w+")
                out_fh.write(msg)
            else:
                print(msg)



    def __init__(self, verbose=False, directed=True, out_degree=3):
        """
        constructor of GBSCA
        :param verbose: if you need the intermediate result
        :param directed: if it's a directed graph
        """
        super().__init__(directed)
        self._VisitedVertex = {}
        self._NotVisitedEdges = set()
        self._cluster_collection = set()
        self._verbose = verbose
        self._cluster_count = 0
        # out degree checking
        self._start_to_checking_outdegree = out_degree

    def set_verbose(self, verbose):
        self._verbose = verbose

    def get_visited_vertex(self):
        # also include their clusters
        return self._VisitedVertex.keys()

    def get_cluster(self, u_element=None):
        if u_element:
            u_vertex = self._vertex_dict[u_element]
            return self._VisitedVertex[u_vertex]
        else:
            return self._cluster_collection

    def get_number_of_clusters(self):
        return self._cluster_count

    def visiting_vertex(self, u, clu):
        self._VisitedVertex[u] = clu
        if u not in clu.get_cluster_nodes().keys():
            raise ValueError("how come")

    def visiting_edge(self, e):
        if isinstance(e, list):
            for i in e:
                if i in self._NotVisitedEdges:
                    self._NotVisitedEdges.remove(i)
        else:
            if e in self._NotVisitedEdges:
                self._NotVisitedEdges.remove(e)

    def insert_edge(self, u, v, **kwargs):
        e = super().insert_edge(u, v, **kwargs)
        if e and e.get_score() > max(kwargs["thresholds"]):
            self._NotVisitedEdges.add(e)

    def remove_low_quality_edges(self, thresholds):
        cutoff = min(thresholds)
        deleted_edges = set()
        for i in self._NotVisitedEdges:
            if i.get_score() <= cutoff:
                deleted_edges.add(i)
                del self._outgoing[i.endpoints()[0]][i.endpoints()[1]]
                del self._incoming[i.endpoints()[1]][i.endpoints()[0]]
        self._NotVisitedEdges -= deleted_edges

    def get_bridge_edges(self, u, clu_v):
        if isinstance(u, self.Vertex):
            return self.get_in_and_out_edges(u, clu_v)

        if isinstance(u, self.Cluster):
            edge_info = {"in_edges": [], "out_edges": [], "double_edges": [],
                         "in_edges_max_score": [], "out_edges_max_score": []}
            for node in u.get_cluster_nodes().keys():
                res = self.get_in_and_out_edges(node, clu_v)
                for i in edge_info.keys():
                    edge_info[i].extend(res[i])
            return edge_info

    def get_in_and_out_edges(self, u, clu_v):
        in_edges = []
        in_edges_max_score = 0
        out_edges = []
        out_edges_max_score = 0
        for x, f in self._outgoing[u].items():
            if x in clu_v.get_cluster_nodes().keys():
                in_edges.append(f)
                in_edges_max_score = max(in_edges_max_score, f.get_score())
        for x, f in self._incoming[u].items():
            if x in clu_v.get_cluster_nodes().keys():
                out_edges.append(f)
                out_edges_max_score = max(out_edges_max_score, f.get_score())

        double_edges = []
        for i in in_edges:
            if self.check_reverse_edge(i):
                double_edges.append(i)
                double_edges.append(self._outgoing[i.endpoints()[1]][i.endpoints()[0]])
        return {"in_edges": in_edges, "out_edges": out_edges, "double_edges": double_edges,
                "in_edges_max_score": [in_edges_max_score if in_edges_max_score != 0 else 1],
                "out_edges_max_score": [out_edges_max_score if out_edges_max_score != 0 else 1]}

    def check_graph_integrity(self, clu):
        key_list = [j for j in self._VisitedVertex.keys() if self._VisitedVertex[j] == clu]
        for i in key_list:
            if i not in clu.get_cluster_nodes().keys():
                raise ValueError("there is a key missing")

        for i in clu.get_cluster_nodes().keys():
            if i not in key_list:
                raise ValueError("there is a key missing")

    def get_candidate_edges(self, edge_threshold):
        print("the size of not visited edges are " + str(len(self._NotVisitedEdges)))
        candidate_edges = [e for e in self._NotVisitedEdges if e.get_score() > edge_threshold]

        return candidate_edges

    def determine_cluster(self, e, cluster_threshold):
        u, v = e.endpoints()
        self.visiting_edge(e)
        method = self.dispatch_determine_method(u, v)
        eval("self."+method+"(u, v, e, cluster_threshold)")

    def dispatch_determine_method(self, u, v):
        cluster_type = [u not in self.get_visited_vertex(), v not in self.get_visited_vertex()]
        if self._verbose:
            msg = [re.sub(r"False", "a cluster", re.sub(r"True", "a node", str(i)))
                   for i in cluster_type]
            print("\t".join(msg))
        if cluster_type == [True, True]:
            return "node_to_node"

        if cluster_type == [True, False]:
            return "node_to_cluster"

        if cluster_type == [False, True]:
            return "cluster_to_node"

        if cluster_type == [False, False]:
            return "cluster_to_cluster"

    def node_to_node(self, u, v, e, cluster_threshold):
        if e.get_score() > cluster_threshold:
            clu = self.Cluster(u, v, e, verbose=self._verbose)
            self._cluster_collection.add(clu)
            self._cluster_count += 1

            self.visiting_vertex(u, clu)
            self.visiting_vertex(v, clu)
        if self._verbose:
            print("successfully initialize a cluster with two seed nodes: " + u.element() + " and " + v.element())

    def node_to_cluster(self, u, v, e, cluster_threshold):
        clu = self._VisitedVertex[v]
        edge_information = self.get_bridge_edges(u, clu)

        if clu.check_and_insert_node(u, cluster_threshold, **edge_information):
            self.visiting_vertex(u, clu)
            self.visiting_edge(edge_information["in_edges"])
            self.visiting_edge(edge_information["out_edges"])
            if self._verbose:
                print("successfully insert the node " + v.element() + " into cluster")
        else:
            if self._verbose:
                print("Unsuccessfully insert the node " + v.element() + " into cluster")

    def cluster_to_node(self, u, v, e, cluster_threshold):
        clu = self._VisitedVertex[u]

        # when the cluster has more nodes than a check point,
        # then only if the the origin of bridge edge, which is in clu cluster,
        # has a outdegree equal to 0, we will process it further
        # the aim of outdegree check is to avoid the situation that two distinct transcripts
        # are falsely grouped together because there is one exons shared between them
        # if clu.get_outdegree(u) != 0 and \
        #        len(clu.get_cluster_nodes().keys()) > self._start_to_checking_outdegree:
        #    return 0

        edge_information = self.get_bridge_edges(v, clu)

        if clu.check_and_insert_node(v, cluster_threshold, **edge_information):
            self.visiting_vertex(v, clu)
            self.visiting_edge(edge_information["in_edges"])
            self.visiting_edge(edge_information["out_edges"])
            if self._verbose:
                print("successfully insert the node " + v.element() + " into cluster")
        else:
            if self._verbose:
                print("Unsuccessfully insert the node " + v.element() + " into cluster")

    def cluster_to_cluster(self, u, v, e, cluster_threshold):
        clu_u = self._VisitedVertex[u]
        clu_v = self._VisitedVertex[v]

        if clu_u == clu_v:
            return 0

        # when the cluster has more nodes than a check point,
        # then only if the the origin of bridge edge, which is in clu_v cluster,
        # has a outdegree equal to 0, we will consider it further
        # if clu_v.get_outdegree(v) != 0 and \
        #        len(clu_v.get_cluster_nodes().keys()) > self._start_to_checking_outdegree:
        #    return 0

        edge_information = self.get_bridge_edges(clu_u, clu_v)

        if clu_u.check_and_merge_cluster(clu_v, cluster_threshold, **edge_information):
            self.visiting_edge(edge_information["in_edges"])
            self.visiting_edge(edge_information["out_edges"])
            for i in clu_v.get_cluster_nodes().keys():
                self._VisitedVertex[i] = clu_u
            self._cluster_collection.remove(clu_v)
            self._cluster_count -= 1

            if self._verbose:
                print("successfully merge two clusters")
        else:
            if self._verbose:
                print("Unsuccessfully merge two clusters")

    def train_part(self, initial_edge_threshold, final_edge_threshold, cluster_threshold):
        edge_threshold = initial_edge_threshold
        lower_amp = (initial_edge_threshold - final_edge_threshold) / 100
        while edge_threshold > final_edge_threshold:
            time_tmp = time.time()
            candidate_edges = list(self.get_candidate_edges(edge_threshold))
            # sort the candidate edges by their scores
            candidate_edges.sort(key=operator.methodcaller('get_score'), reverse=True)
            print("first part takes " + str(time.time() - time_tmp))
            time_tmp = time.time()
            while candidate_edges:
                # process the edges one by one
                test_edge = candidate_edges.pop(0)
                self.determine_cluster(test_edge, cluster_threshold)

            print("second part takes " + str(time.time() - time_tmp))
            edge_threshold = edge_threshold - lower_amp
            print(self.get_number_of_clusters())
            print("edge_threshold = " + str(edge_threshold))
            print()

    def print_cluster_report(self, out_base=None):
        cnt = 0
        clu_nodes = []
        output_cluster = open(out_base + "_clusters.txt", "w")
        for clu in self._cluster_collection:
            cnt += 1
            clu_nodes.append(len(clu.get_cluster_nodes().keys()))

            output_cluster.write("\n".join([str(i) + "\tCluster_" + str(cnt)
                                            for i in
                                            map(operator.methodcaller('element'), clu.get_cluster_nodes().keys())])
                                 + "\n")

        print("there are " + str(cnt) + " clusters.\nTheir number of nodes are "
              + ", ".join([str(i) for i in clu_nodes]) +
              "\n" + str(sum(clu_nodes)) + " nodes are involved in the clusters")

    def check_with_known_isoforms(self):
        luci_gene_id = "TRINITY_DN8457_c0_g1"
        ribo_gene_id = "TRINITY_DN3636_c0_g1"
        check_luci = set()
        check_ribo = set()
        for clu in self._cluster_collection:
            for i in map(operator.methodcaller('element'), clu.get_cluster_nodes().keys()):
                if re.match(luci_gene_id, i):
                    check_luci.add(clu)
                if re.match(ribo_gene_id, i):
                    check_ribo.add(clu)

        print("the luciferases are in " + str(len(check_luci)) + " different clusters")
        print("the ribosomals are in " + str(len(check_ribo)) + " different clusters")
        return check_luci, check_ribo

    def plot_a_cluster(self, clu):
        nodes_set = clu.get_cluster_nodes().keys()
        G = networkx.DiGraph()
        G.add_nodes_from(nodes_set)

        edges_set = []
        for i in nodes_set:
            edges_set.extend([e.endpoints() for e in self.incident_edges(i) if e.endpoints()[1] in nodes_set])
        G.add_edges_from(edges_set)
        networkx.draw_circular(G, labels={i: re.sub(r"TRINITY_DN", "", i.element()) for i in nodes_set})

    def plot_adjacent(self, node):
        node_obj = self._vertex_dict[node]

        G = networkx.DiGraph()
        incident_node_set = set()
        for i in self.incident_vertex(node_obj, outgoing=True):
            incident_node_set.add(i)

        for i in self.incident_vertex(node_obj, outgoing=False):
            incident_node_set.add(i)

        G.add_nodes_from(incident_node_set)

        edges_col = {}
        for e in self._outgoing[node_obj].values():
            G.add_edge(e.endpoints()[0], e.endpoints()[1], weights=round(e.get_score(), 3), color='r')
            edges_col[e] = "r"

        for e in self._incoming[node_obj].values():
            G.add_edge(e.endpoints()[0], e.endpoints()[1], weights=round(e.get_score(), 3), color='b')
            edges_col[e] = "b"

        pos = networkx.circular_layout(G)  # positions for all nodes

        # nodes
        networkx.draw_networkx_nodes(G, pos, node_size=200)

        # edges
        networkx.draw_networkx_edges(G, pos, width=2)
        networkx.draw_networkx_edge_labels(G, pos, edge_labels=networkx.get_edge_attributes(G, 'weights'))

        # labels
        networkx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif',
                                      labels={i: re.sub(r"TRINITY_DN", "", i.element()) for i in incident_node_set})


def parse_blast_outfmt(outfmt):
    """
    :param outfmt: the output format specified during blast
    :return: the location of column qseqid, sseqid, pident, qcovs
    """
    column = re.split(r"\s+", outfmt)
    format_number = [i for i in column if re.match(r"^\d+$", i)]
    if format_number == ["6"] or format_number == []:
        init = 1 if format_number == ["6"] else 0
        try:
            qseqid = column.index('qseqid') - init
            sseqid = column.index('sseqid') - init
            pident = column.index('pident') - init
            qcovs = column.index('qcovs') - init
            evalue = column.index('evalue') - init
            return [qseqid, sseqid, pident, qcovs, evalue]
        except ValueError as err:
            print('please specify valid blast column id. the error message is: ', err)
    else:
        raise ValueError("outfmt should be format 6 with appropriate column names")


def main(blast_id, blast_dir, out_base, outfmt, use_evalue,
         cluster_threshold, initial_edge_threshold, final_edge_threshold):

    start_time = time.time()
    qseqid, sseqid, pident, qcovs, evalue = parse_blast_outfmt(outfmt)
    if not use_evalue:
        evalue = None

    # start to construct the graph
    my_graph = GBSCA()
    my_graph.read_vertex(blast_id)
    my_graph.read_edge(blast_dir, qseqid, sseqid, pident, qcovs, evalue, [cluster_threshold, final_edge_threshold])
    print("Graph initialization is finished.")

    # begin with the most stringent threshold and gradually lower it down until it hit the final threshold
    my_graph.train_part(initial_edge_threshold, final_edge_threshold, cluster_threshold)

    my_graph.print_cluster_report(out_base)
    my_graph.check_with_known_isoforms()
    print("Overall time is " + str(time.time() - start_time))

    with open(out_base + "_GBSCA.obj", 'wb') as clusters_file:
        pickle.dump(my_graph, clusters_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--blastID', required=True,
                        help='The input file directory where your blast id are stored')
    parser.add_argument('--blastDir', required=True,
                        help='The input file directory where your blast result are stored')
    parser.add_argument('--outBase', required=True,
                        help='The base name for the output files')
    parser.add_argument('--outfmt',
                        default='6 qseqid sseqid bitscore qcovs evalue pident qstart qend sstart send length qlen',
                        help='blast output format')
    parser.add_argument('--use_evalue', action="store_true", help='')
    parser.add_argument('--cluster_threshold', default=0.6, type=float,
                        help='used to determine whether a node can be added to a cluster or two clusters can be merged')
    parser.add_argument('--initial_edge_threshold', default=0.98, type=float,
                        help='filter candidate edges that serve as seed')
    parser.add_argument('--final_edge_threshold', default=0.8, type=float,
                        help='filter candidate edges that serve as seed')
    parser.add_argument('--verbose',
                        action="store_true",
                        help='if we need to print out some intermediate result')

    argument = parser.parse_args()
    main(argument.blastID, argument.blastDir, argument.outBase, argument.outfmt, argument.use_evalue,
         argument.cluster_threshold, argument.initial_edge_threshold, argument.final_edge_threshold)






