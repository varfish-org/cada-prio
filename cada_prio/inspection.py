"""Helpers for inspection models"""

import pickle

import networkx as nx

from cada_prio.predict import load_hgnc_info


def dump_graph(path_graph: str, path_hgnc_info: str):
    _, hgnc_info_by_id = load_hgnc_info(path_hgnc_info)
    with open(path_graph, "rb") as inputf:
        graph: nx.Graph = pickle.load(inputf)
    for edge in sorted(graph.edges):
        lhs, rhs = edge
        if lhs.startswith("HGNC:"):
            lhs = "Entrez:%s" % hgnc_info_by_id[lhs].ncbi_gene_id
        if rhs.startswith("HGNC:"):
            rhs = "Entrez:%s" % hgnc_info_by_id[rhs].ncbi_gene_id
        print(f"{lhs}\t{rhs}")
