"""carmm vib backend """

from pathlib import Path

import numpy as np

# from click 7.1 on
_default_context_settings = {"show_default": True}

@click.command(cls=ClickAliasedGroup)
def vib():
    """perform vib analysis"""

@vib.command()
@click.argument("file")
@click.option("-p", "--plot", default=False, help="Plot vibrational analysis")

def vib_analysis(file):
import carmm.analyse.vibrations as vib



