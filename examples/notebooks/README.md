[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/logsdail/carmm/HEAD?labpath=examples%2Fnotebooks)

This is folder with Jupyter notebooks explaining how to use CARMM packages.

If you want to play with the notebooks, run on your local computer with `jupyter notebook` and then selecting the notebook of interest.

(*Note: you might need to install jupyter first with `pip install --user jupyter`*)

Be sure to set the environment paths appropriately so you can load modules:

On Linux: set `JUPYTER_PATH` variable in a similar way that `PYTHONPATH` is defined [previously](../../README.md).

**Example:** `export JUPYTER_PATH=/path/to/carmm/folder:$JUPYTER_PATH`

On Windows (also works on Linux): Add the appropriate path to the environment in the first cell of the notebook

**Example:** `import sys; sys.path.append(r'DriveLetter:\path\to\carmm\')`


Alternatively, clicking on the binder link above will take you to a web platform where you can play with the notebooks without installing anything!
