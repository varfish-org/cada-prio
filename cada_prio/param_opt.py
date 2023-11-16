"""Code for parameter optimization."""

from datetime import timedelta
import gzip
import hashlib
import json
import os
import random
import time
import typing

import cattrs
from logzero import logger
import tqdm

from cada_prio import predict, train_model

#: Default seeed value.
DEFAULT_SEED = 1


def load_clinvar_phenotype_links_jsonl(path) -> typing.Iterable[train_model.GenePhenotypeRecord]:
    if path.endswith(".gz"):
        inputf = gzip.open(path, "rt")
    else:
        inputf = open(path, "rt")
    with inputf:
        for line in inputf:
            yield cattrs.structure(json.loads(line), train_model.GenePhenotypeRecord)


def make_links_unique_by_submitter(
    links: typing.Iterable[train_model.GenePhenotypeRecord],
) -> typing.List[train_model.GenePhenotypeRecord]:
    links_dict = {
        f"{link.submitter}/{','.join(link.hgnc_ids)}/{','.join(link.hpo_terms)}": link
        for link in links
        if len(link.hgnc_ids) == 1
    }
    return list(links_dict.values())


def load_links(
    path_clinvar_phenotype_links: str,
    fraction_links: typing.Optional[float],
    path_validation_links: typing.Optional[str],
    seed: typing.Optional[int],
):
    rng = random.Random(seed or DEFAULT_SEED)

    logger.info("Loading phenotype links...")
    logger.info("- load JSONL")
    orig_phenotype_links = list(load_clinvar_phenotype_links_jsonl(path_clinvar_phenotype_links))
    logger.info("- make unique by submitter")
    phenotype_links = make_links_unique_by_submitter(orig_phenotype_links)
    logger.info("- randomizing")
    rng.shuffle(phenotype_links)
    logger.info("... done loading %s links", len(phenotype_links))

    if fraction_links is None:
        assert (
            path_validation_links is not None
        ), "must give fraction_links or path_validation_links"
        links_training = phenotype_links
        orig_links_validation = list(load_clinvar_phenotype_links_jsonl(path_validation_links))
        links_validation = make_links_unique_by_submitter(orig_links_validation)
        links_test: typing.List[train_model.GenePhenotypeRecord] = []
        logger.info("Counts in explicit validation set...")
        logger.info("- training:   % 6d", len(links_training))
        logger.info(
            "- validation: % 6d (non-unique: %d)", len(links_validation), len(orig_links_validation)
        )
        logger.info("- test:       % 6d", len(links_test))
        logger.info("... that's all")
    else:
        logger.info("Computing counts...")
        n_links_total = len(phenotype_links)
        n_links_training_all = int(0.6 * n_links_total)
        n_links_training = int(n_links_training_all * fraction_links)
        n_links_validation = (n_links_total - n_links_training_all) // 2
        n_links_test = n_links_validation
        links_training = phenotype_links[:n_links_training]
        links_validation = phenotype_links[n_links_training : n_links_training + n_links_validation]
        links_test = phenotype_links[
            n_links_training
            + n_links_validation : n_links_training
            + n_links_validation
            + n_links_test
        ]
        logger.info("- total:        % 6d", len(phenotype_links))
        logger.info("- training-all: % 6d", n_links_training_all)
        logger.info("- training:     % 6d", len(links_training))
        logger.info("- validation:   % 6d", len(links_validation))
        logger.info("- test:         % 6d", len(links_test))
        logger.info("... done computing counts")

    return links_training, links_validation, links_test


def prepare_training(
    path_out: str,
    links_training: typing.List[train_model.GenePhenotypeRecord],
    links_validation: typing.List[train_model.GenePhenotypeRecord],
    links_test: typing.List[train_model.GenePhenotypeRecord],
    fraction_links: typing.Optional[float],
    path_validation_links: typing.Optional[str],
    path_embedding_params: typing.Optional[str],
    seed: typing.Optional[int],
):
    if seed is None:
        seed = DEFAULT_SEED

    logger.info("Preparing training...")
    if path_embedding_params:
        logger.info("- loading embedding params from %s", path_embedding_params)
        with open(path_embedding_params, "rt") as inputf:
            embedding_params = cattrs.structure(json.load(inputf), train_model.EmbeddingParams)
    else:
        logger.info("- using default embedding params params")
        embedding_params = train_model.EmbeddingParams(seed_embedding=seed + 23, seed_fit=seed + 42)
    params_name = hashlib.md5(
        json.dumps(cattrs.unstructure(embedding_params)).encode("utf-8")
    ).hexdigest()

    if fraction_links is None:
        path_out = f"{path_out}/explicit-validation.{seed}.{params_name}"
        _ = path_validation_links
    else:
        path_out = f"{path_out}/{int(fraction_links * 100)}.{seed}.{params_name}"

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
    path_embedding_params: typing.Optional[str],
    cpus: int,
):
    start = time.time()
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
    elapsed = time.time() - start
    logger.info("... done running training in %s", timedelta(seconds=elapsed))


def run_validation(
    path_out: str,
    path_links_validation: str,
) -> typing.Dict[int, float]:
    """Run validation on the links from the path.

    :return: dictionary with top N -> percentage of links in top N
    """
    start = time.time()
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

    result = {
        no: 100.0 * count / links_total if links_total > 0 else 0.0
        for no, count in links_top.items()
    }
    logger.info("<result>")
    print(json.dumps(result, indent=2))
    logger.info("</result>")

    elapsed = time.time() - start
    logger.info("... done running validation in %s", timedelta(seconds=elapsed))
    return result


def train_and_validate(
    *,
    path_out: str,
    path_hgnc_json: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_clinvar_phenotype_links: str,
    fraction_links: typing.Optional[float],
    path_validation_links: typing.Optional[str],
    path_embedding_params: typing.Optional[str],
    seed: typing.Optional[int],
    cpus: int,
) -> typing.Dict[int, float]:
    """Train model and run validation on the links from the path.

    :return: dictionary with top N -> percentage of links in top N
    """

    links_training, links_validation, links_test = load_links(
        path_clinvar_phenotype_links, fraction_links, path_validation_links, seed
    )

    (
        path_out,
        path_links_training,
        path_links_validation,
        path_links_test,
        path_embedding_params,
    ) = prepare_training(
        path_out,
        links_training,
        links_validation,
        links_test,
        fraction_links,
        path_validation_links,
        path_embedding_params,
        seed,
    )
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

    return run_validation(
        path_out,
        path_links_validation,
    )
