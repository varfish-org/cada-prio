"""Code for parameter optimization."""

import gzip
import hashlib
import json
import os
import random
import typing

import cattrs
from logzero import logger
import tqdm

from cada_prio import predict, train_model


def load_clinvar_phenotype_links_jsonl(path) -> typing.Iterable[train_model.GenePhenotypeRecord]:
    if path.endswith(".gz"):
        inputf = gzip.open(path, "rt")
    else:
        inputf = open(path, "rt")
    with inputf:
        for line in inputf:
            yield cattrs.structure(json.loads(line), train_model.GenePhenotypeRecord)


def load_links(path_clinvar_phenotype_links: str, fraction_links: float, seed: int):
    rng = random.Random(seed)

    logger.info("Loading phenotype links...")
    logger.info("- load JSONL")
    orig_phenotype_links = list(load_clinvar_phenotype_links_jsonl(path_clinvar_phenotype_links))
    logger.info("- make unique by submitter")
    phenotype_links_dict = {
        f"{link.submitter}/{','.join(link.hgnc_ids)}/{','.join(link.hpo_terms)}": link
        for link in orig_phenotype_links
        if len(link.hgnc_ids) == 1
    }
    phenotype_links = list(phenotype_links_dict.values())
    logger.info("- randomizing")
    rng.shuffle(phenotype_links)
    logger.info("... done loading %s links", len(phenotype_links))

    logger.info("Computing counts...")
    n_links_total = len(phenotype_links)
    n_links_used = int(fraction_links * n_links_total)
    n_links_training = int(n_links_used * 0.6)
    n_links_validation = (n_links_used - n_links_training) // 2
    n_links_test = n_links_validation
    links_training = phenotype_links[:n_links_training]
    links_validation = phenotype_links[n_links_training : n_links_training + n_links_validation]
    links_test = phenotype_links[
        n_links_training + n_links_validation : n_links_training + n_links_validation + n_links_test
    ]
    logger.info("- total:      % 6d", len(phenotype_links))
    logger.info("- used:       % 6d", n_links_used)
    logger.info("- training:   % 6d", len(links_training))
    logger.info("- validation: % 6d", len(links_validation))
    logger.info("- test:       % 6d", len(links_test))
    logger.info("... done computing counts")

    return links_training, links_validation, links_test


def prepare_training(
    links_training: typing.List[train_model.GenePhenotypeRecord],
    links_validation: typing.List[train_model.GenePhenotypeRecord],
    links_test: typing.List[train_model.GenePhenotypeRecord],
    fraction_links: float,
    seed: int,
):
    logger.info("Preparing training...")
    embedding_params = train_model.EmbeddingParams(seed_embedding=seed + 23, seed_fit=seed + 42)
    params_name = hashlib.md5(
        json.dumps(cattrs.unstructure(embedding_params)).encode("utf-8")
    ).hexdigest()
    path_out = f"param_opt.d/{int(fraction_links * 100)}.{seed}.{params_name}"
    os.makedirs(path_out, exist_ok=True)
    logger.info("- output path: %s", path_out)
    logger.info(
        "- embedding params:\n%s", json.dumps(cattrs.unstructure(embedding_params), indent=2)
    )
    path_embedding_params = f"{path_out}/embedding_params.json"
    with open(path_embedding_params, "wt") as outputf:
        print(json.dumps(cattrs.unstructure(embedding_params), indent=2), file=outputf)
    path_links_training = f"{path_out}/links_training.jsonl"
    logger.info("- trailing links: %s", path_links_training)
    with open(path_links_training, "wt") as outputf:
        for link in links_training:
            print(json.dumps(cattrs.unstructure(link)), file=outputf)
    path_links_validation = f"{path_out}/links_validation.jsonl"
    logger.info("- validation links: %s", path_links_validation)
    with open(path_links_validation, "wt") as outputf:
        for link in links_validation:
            print(json.dumps(cattrs.unstructure(link)), file=outputf)
    path_links_test = f"{path_out}/links_test.jsonl"
    logger.info("- test links: %s", path_links_test)
    with open(path_links_test, "wt") as outputf:
        for link in links_test:
            print(json.dumps(cattrs.unstructure(link)), file=outputf)
    logger.info("... done preparing training")

    return (
        path_out,
        path_links_training,
        path_links_validation,
        path_links_test,
        path_embedding_params,
    )


def run_training(
    path_out: str,
    path_hgnc_json: str,
    path_links_training: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_embedding_params: str,
    cpus: int,
):
    logger.info("Running training...")
    train_model.run(
        f"{path_out}/model",
        path_hgnc_json,
        path_links_training,
        path_hpo_genes_to_phenotype,
        path_hpo_obo,
        path_embedding_params,
        cpus=cpus,
    )
    logger.info("... done running training")


def run_validation(
    path_out: str,
    path_links_validation: str,
):
    logger.info("Running validation...")
    path_model = f"{path_out}/model"
    logger.info("- model: %s", path_model)
    all_to_hgnc, _ = predict.load_hgnc_info(os.path.join(path_model, "hgnc_info.jsonl"))
    graph, model = predict.load_graph_model(path_model)
    _, hpo_id_from_alt, _ = train_model.load_hpo_ontology(os.path.join(path_model, "hp.obo"))

    logger.info("... validation steps ...")
    links_total = 0
    links_top = {
        1: 0,
        5: 0,
        10: 0,
        50: 0,
        100: 0,
    }
    pb = tqdm.tqdm(list(load_clinvar_phenotype_links_jsonl(path_links_validation)), unit=" records")
    for link_validation in pb:
        if len(link_validation.hgnc_ids) != 1:
            logger.warn(
                "skipping submission %s with %d genes",
                link_validation.scv,
                len(link_validation.hgnc_ids),
            )
            continue
        else:
            hgnc_id = link_validation.hgnc_ids[0]

        try:
            hpo_terms = list(
                sorted(set(hpo_id_from_alt.get(x, x) for x in link_validation.hpo_terms))
            )
            _, sorted_scores = predict.run_prediction(hpo_terms, None, all_to_hgnc, graph, model)
        except predict.NoValidHpoTerms:
            logger.warn("no valid HPO terms in %s (skipped)", link_validation.hpo_terms)
            continue

        links_total += 1
        ranked_genes = list(sorted_scores.keys())
        if hgnc_id in ranked_genes:
            rank = ranked_genes.index(hgnc_id) + 1
        else:
            rank = len(ranked_genes)
        for no in links_top.keys():
            if rank <= no:
                links_top[no] += 1

    result = {no: 100.0 * count / links_total for no, count in links_top.items()}
    logger.info("<result>")
    print(json.dumps(result, indent=2))
    logger.info("</result>")

    logger.info("... done running validation")


def perform_parameter_optimization(
    path_hgnc_json: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_clinvar_phenotype_links: str,
    fraction_links: float,
    seed: int,
    cpus: int,
):
    """Simulate cases based on the dataset file."""

    links_training, links_validation, links_test = load_links(
        path_clinvar_phenotype_links, fraction_links, seed
    )

    (
        path_out,
        path_links_training,
        path_links_validation,
        path_links_test,
        path_embedding_params,
    ) = prepare_training(links_training, links_validation, links_test, fraction_links, seed)
    _ = path_links_test

    run_training(
        path_out,
        path_hgnc_json,
        path_links_training,
        path_hpo_genes_to_phenotype,
        path_hpo_obo,
        path_embedding_params,
        cpus=cpus,
    )

    run_validation(
        path_out,
        path_links_validation,
    )
