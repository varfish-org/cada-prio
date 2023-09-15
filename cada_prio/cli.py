"""Console script for CADA"""

import logging
import os
import sys
import typing

import click
import logzero

from cada_prio import _version, inspection, param_opt, predict, train_model

# Lower the update interval of tqdm to 5 seconds if stdout is not a TTY.
if not sys.stdout.isatty():
    os.environ["TQDM_MININTERVAL"] = "5"


@click.group()
@click.version_option(_version.__version__)
@click.option("--verbose/--no-verbose", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    if verbose:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)


@cli.command("train")
@click.argument("path_out", type=str)
@click.option("--path-hgnc-json", type=str, help="path to HGNC JSON", required=True)
@click.option(
    "--path-gene-hpo-links", type=str, help="path to gene-HPO term links JSONL", required=True
)
@click.option(
    "--path-hpo-genes-to-phenotype",
    type=str,
    help="path to genes_to_phenotype.txt file",
    required=True,
)
@click.option("--path-hpo-obo", type=str, help="path HPO OBO file", required=True)
@click.option(
    "--path-embedding-params", type=str, help="optional path to JSON file with embedding parameters"
)
@click.option("--cpus", type=int, help="number of CPUs to use", default=1)
@click.pass_context
def cli_train(
    ctx: click.Context,
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_embedding_params: typing.Optional[str],
    cpus: int,
):
    """train model"""
    ctx.ensure_object(dict)
    train_model.run(
        path_out,
        path_hgnc_json,
        path_gene_hpo_links,
        path_hpo_genes_to_phenotype,
        path_hpo_obo,
        path_embedding_params,
        cpus,
    )


@cli.command("predict")
@click.argument("path_model", type=str)
@click.option(
    "--hpo-terms",
    type=str,
    help="comma-separate HPO terms or @file with space-separated ones",
    required=True,
)
@click.option(
    "--genes",
    type=str,
    help="comma-separated genes to restrict prediction to or @file with space-separated ones",
)
@click.pass_context
def cli_predict(
    ctx: click.Context,
    path_model: str,
    hpo_terms: str,
    genes: typing.Optional[str],
):
    """perform prediction/prioritication"""
    ctx.ensure_object(dict)

    if hpo_terms.startswith("@"):
        with open(hpo_terms[1:]) as f:
            hpo_term_list = [x.strip() for x in f.read().split()]
    else:
        hpo_term_list = [x.strip() for x in hpo_terms.split(",")]

    gene_list = None
    if genes:
        if genes.startswith("@"):
            with open(genes[1:]) as f:
                gene_list = [x.strip() for x in f.read().split()]
        else:
            gene_list = [x.strip() for x in genes.split(",")]

    if predict.run(path_model, hpo_term_list, gene_list) != 0:
        ctx.exit(1)


@cli.group("utils")
def cli_utils():
    """utilities"""


@cli_utils.command("dump-graph")
@click.argument("path_graph", type=str)
@click.argument("path_hgnc_info", type=str)
@click.pass_context
def cli_dump_graph(
    ctx: click.Context,
    path_graph: str,
    path_hgnc_info: str,
):
    """dump graph edges for debugging"""
    ctx.ensure_object(dict)
    inspection.dump_graph(path_graph, path_hgnc_info)


@cli.group("tune")
def cli_tune():
    """hyperparameter tuning"""


@cli_tune.command("train-eval")
@click.argument("path_out", type=str)
@click.option("--path-hgnc-json", type=str, help="path to HGNC JSON", required=True)
@click.option(
    "--path-hpo-genes-to-phenotype",
    type=str,
    help="path to genes_to_phenotype.txt file",
    required=True,
)
@click.option("--path-hpo-obo", type=str, help="path HPO OBO file", required=True)
@click.option(
    "--path-clinvar-phenotype-links",
    type=str,
    help="path to ClinVar phenotype links JSONL",
    required=True,
)
@click.option(
    "--fraction-links",
    type=float,
    help="fraction of links to add to the graph (conflicts with --path-validation-links)",
)
@click.option(
    "--path-validation-links",
    type=str,
    help="path to validation links JSONL (conflicts with --fraction-links)",
)
@click.option(
    "--path-embedding-params", type=str, help="path to JSON file with embedding params; optional"
)
@click.option(
    "--seed",
    type=int,
    help="seed for random number generator",
)
@click.option("--cpus", type=int, help="number of CPUs to use", default=1)
@click.pass_context
def cli_param_opt(
    ctx: click.Context,
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
):
    """train and evaluate model for one set of parameters"""
    if bool(fraction_links) == bool(path_validation_links):
        raise click.UsageError(
            "exactly one of --fraction-links and --path-validation-links must be given"
        )

    ctx.ensure_object(dict)
    param_opt.train_and_validate(
        path_out=path_out,
        path_hgnc_json=path_hgnc_json,
        path_hpo_genes_to_phenotype=path_hpo_genes_to_phenotype,
        path_hpo_obo=path_hpo_obo,
        path_clinvar_phenotype_links=path_clinvar_phenotype_links,
        fraction_links=fraction_links,
        path_validation_links=path_validation_links,
        path_embedding_params=path_embedding_params,
        seed=seed,
        cpus=cpus,
    )
