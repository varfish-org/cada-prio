"""Model-based prediction"""

import json
import os
import typing

import cattrs
from gensim.models import Word2Vec
from logzero import logger
import networkx as nx
import numpy as np

from cada_prio import train_model


def load_hgnc_info(path_hgnc_json: str) -> typing.List[train_model.GeneIds]:
    result = []
    with open(path_hgnc_json, "rt") as f:
        for line in f:
            result.append(cattrs.structure(json.loads(line), train_model.GeneIds))
    return result


def run(
    path_model: str,
    orig_hpo_terms: typing.List[str],
    genes: typing.Optional[typing.List[str]] = None,
) -> int:
    # Load and prepare data
    logger.info("Loading HGNC info...")
    logger.info("- parsing")
    hgnc_infos = load_hgnc_info(os.path.join(path_model, "hgnc_info.jsonl"))
    logger.info("- create mapping")
    all_to_hgnc = {}
    for record in hgnc_infos:
        all_to_hgnc[record.symbol] = record
        all_to_hgnc[record.ncbi_gene_id] = record
        all_to_hgnc[record.hgnc_id] = record
        if record.ensembl_gene_id:
            all_to_hgnc[record.ensembl_gene_id] = record
    hgnc_info_by_id = {record.hgnc_id: record for record in hgnc_infos}
    hgnc_ids = []
    for gene in genes or []:
        if gene in all_to_hgnc:
            hgnc_ids.append(all_to_hgnc[gene].hgnc_id)
        else:
            logger.warning("could not resolve HGNC ID for gene %s", gene)
    logger.info("... done loading HGNC info")

    logger.info("Loading graph...")
    graph = nx.read_gpickle(os.path.join(path_model, "graph.gpickle"))
    logger.info("... done loading graph")
    logger.info("Loading model...")
    model = Word2Vec.load(os.path.join(path_model, "model"))
    logger.info("... done loading model")

    # Lookup HPO term embeddings.
    hpo_terms = {}
    for hpo_term in orig_hpo_terms:
        if hpo_term not in model.wv:
            logger.warn("skipping HPO term %s as it is not in the model", hpo_term)
        else:
            hpo_terms[hpo_term] = model.wv[hpo_term]
    if not hpo_terms:
        logger.error("no valid HPO terms in model")
        return 1

    # Generate a score for each gene in the knowledge graph
    logger.info("Generating scores...")
    gene_scores = {}
    for node_id in graph.nodes():
        if node_id.startswith("HGNC:"):  # is gene
            hgnc_id = node_id
            if genes and hgnc_id not in genes:
                continue  # skip, not asked for

            this_gene_scores = []
            hgnc_id_emb = model.wv[hgnc_id]
            for hpo_term, hpo_term_emb in hpo_terms.items():
                score = np.dot(hpo_term_emb, hgnc_id_emb)
                this_gene_scores.append(score)
            gene_scores[hgnc_id] = sum(this_gene_scores) / len(hpo_terms)

    # Write out results to stdout, largest score first
    sorted_scores = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
    print("# query (len=%d): %s" % (len(hpo_terms), ",".join(hpo_terms)))
    print("\t".join(["rank", "score", "gene_symbol", "ncbi_gene_id", "hgnc_id"]))
    for rank, (hgnc_id, score) in enumerate(sorted_scores, start=1):
        hgnc_info = hgnc_info_by_id[hgnc_id]
        print(
            "\t".join(
                map(
                    str,
                    [
                        rank,
                        score,
                        hgnc_info.symbol,
                        hgnc_info.ncbi_gene_id,
                        hgnc_info.hgnc_id,
                    ],
                )
            )
        )

    return 0
