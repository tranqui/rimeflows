## Overview

This package contains the calculations underlying the following work:

* J. F. Robinson, P. B. Warren, M. R. Turner and R. P. Sear, "Critical inertia for particle capture is determined by surface geometry at forward stagnation point", [arXiv: 2310.03474](https://arxiv.org/abs/2310.03474).

This code performs calculations of point particle capture efficiency in model 2d flows. This data is used to render the figures in the manuscript above. Specifically, the code underlying this manuscript is found in (number corresponds to each figure):

1. [`sketch-ellipse.py`](https://github.com/tranqui/rimeflows/blob/main/sketch-ellipse.py)
2. [`show-kuwabara.py`](https://github.com/tranqui/rimeflows/blob/main/show-kuwabara.py)
3. [`phase-portraits.py`](https://github.com/tranqui/rimeflows/blob/main/phase-portraits.py)
4. [`show-singular-onset.py`](https://github.com/tranqui/rimeflows/blob/main/show-singular-onset.py)
5. [`unknown.ipynb`](https://github.com/tranqui/rimeflows/blob/main/unknown.ipynb)
6. [`show-singular-onset.py`](https://github.com/tranqui/rimeflows/blob/main/show-singular-onset.py)
7. [`unknown.ipynb`](https://github.com/tranqui/rimeflows/blob/main/unknown.ipynb)

These scripts primarily plot data that was pregenerated using the scripts

1. [`critical-stokes.py`](https://github.com/tranqui/rimeflows/blob/main/critical-stokes.py): measures the critical Stokes number above which particles are deposited onto a surface.
2. [`efficiency.py`](https://github.com/tranqui/rimeflows/blob/main/efficiency.py): measures the efficiency of particle deposition above the critical Stokes number.

These two scripts work with the flow fields which are defined in the [`flows`](https://github.com/tranqui/rimeflows/blob/main/flows.py) subdirectory. The precalculated data used in the scripts was obtained via executing the two scripts above with various argument through the script [`batch-calcs.sh`](https://github.com/tranqui/rimeflows/blob/main/batch-calcs.sh).

### Copying

These programs are free software: you can redistribute them and/or modify
them under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with these programs.  If not, see
<http://www.gnu.org/licenses/>.


### Copyright

These programs are copyright &copy; 2025 Joshua Robinson and copyright
&copy; 2025 Patrick B Warren and copyright
&copy; 2025 Richard Sear, as indicated in the individual files.
