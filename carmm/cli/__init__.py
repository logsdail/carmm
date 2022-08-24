"""CLI for carmm with click"""

import shutil
import click
import click_completion


from . import analyse
click_completion.init()


@click.command(cls=AliasedGroup)
@click.version_option(vibes_version, "-V", "--version")
@click.option("-v", "--verbose", is_flag=True, hidden=True)
@click.option("--silent", is_flag=True, help="set verbosity level to 0")
@click.pass_context

def cli(ctx, verbose, silent):
    "carmm: "
    if verbose:
        verbosity = 2
    elif silent:
        verbosity = 0
    else:
        verbosity = 1

    ctx.obj = CliTracker(verbose=verbosity)
    ctx.help_option_names = ["-h", "--help"]

    if verbosity > 1:
        click.echo(f"Welcome to carmm!\n")


cli.add_command(analyse.analyse)



