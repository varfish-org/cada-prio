"""Console script for CADA"""

import typing

import click

from cada_prio import _version, predict, train_model


@click.group()
@click.version_option(_version.__version__)
@click.option("--verbose/--no-verbose", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command("train-model")
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
@click.option("--cpus", type=int, help="number of CPUs to use", default=1)
@click.pass_context
def cli_train_model(
    ctx: click.Context,
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
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
