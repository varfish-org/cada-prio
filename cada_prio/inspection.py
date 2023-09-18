"""Helpers for inspection models"""

import pickle

from logzero import logger
import networkx as nx

from cada_prio.predict import load_hgnc_info


def dump_graph(path_graph: str, path_hgnc_info: str, hgnc_to_entrez: bool):
    _, hgnc_info_by_id = load_hgnc_info(path_hgnc_info)
    hgnc_info_by_ncbi_gene_id = {
        "Entrez:%s" % x.ncbi_gene_id: x for x in hgnc_info_by_id.values() if x.ncbi_gene_id
    }
    with open(path_graph, "rb") as inputf:
        graph: nx.Graph = pickle.load(inputf)
    for edge in sorted(graph.edges):
        lhs, rhs = edge

        if hgnc_to_entrez:
            if lhs.startswith("HGNC:"):
                lhs = "Entrez:%s" % hgnc_info_by_id[lhs].ncbi_gene_id
            if rhs.startswith("HGNC:"):
                rhs = "Entrez:%s" % hgnc_info_by_id[rhs].ncbi_gene_id
        else:  # entrez to hgnc
            if lhs.startswith("Entrez:"):
                if lhs not in hgnc_info_by_ncbi_gene_id:
                    logger.warning("no HGNC ID for Entrez ID %s", lhs)
                    continue
                lhs = hgnc_info_by_ncbi_gene_id[lhs].hgnc_id
            if rhs.startswith("Entrez:"):
                if rhs not in hgnc_info_by_ncbi_gene_id:
                    logger.warning("no HGNC ID for Entrez ID %s", rhs)
                    continue
                rhs = hgnc_info_by_ncbi_gene_id[rhs].hgnc_id

        print(f"{lhs}\t{rhs}")
