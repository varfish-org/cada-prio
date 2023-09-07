"""Console script for CADA"""

import click

from cada_prio import _version


@click.group()
@click.version_option(_version.__version__)
@click.option("--verbose/--no-verbose", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
