[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/logsdail/carmm/master?filepath=examples%2Fnotebooks)

This is folder with Jupyter notebooks explaining how to use CARMM packages.

If you want to play with the notebooks, run on your local computer with `jupyter notebook` and then selecting the notebook of interest.

(*Note: you might need to install jupyter first with `pip install --user jupyter`*)

Be sure to set the `JUPYTER_PATH` variable in a similar way that `PYTHONPATH` is defined [previously](../../README.md).

*Example: export JUPYTER_PATH=/path/to/carmm/folder:$JUPYTER_PATH*

Alternatively, clicking on the binder link above will take you to a web platform where you can play with the notebooks without installing anything!

In case you are still struggling to access and run this notebook, please do the following at the start of the notebook. 
```
import sys
sys.path.append(r'C:\path\to\your\carmm\folder\on\your\system')
import carmm
```
