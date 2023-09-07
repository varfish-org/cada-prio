"""Console script for CADA"""

import click

from cada_prio import _version, train_model


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
@click.pass_context
def cli_train_model(
    ctx: click.Context,
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
):
    ctx.ensure_object(dict)
    train_model.run(
        path_out, path_hgnc_json, path_gene_hpo_links, path_hpo_genes_to_phenotype, path_hpo_obo
    )
