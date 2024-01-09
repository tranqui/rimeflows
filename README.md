## Installation

Prerequisites:
* ?

## Calculations of point particle capture efficiency in model 2d flows


### Main scripts

There are several flow fields implemented in this package. We provide the script `show-flow.py` to visualise the various flow fields. The Kuwabara flow field is the most important for flow through porous media. Visualise this flow field by running e.g. (the parameters simply control the resulting colors):

    python show-flow.py 0.1 -cc grey -fc '#e6e6ff' -sc steelblue -ctc red -ntc blue

Which should generate this figure:

![Kuwabara particle trajectories](https://github.com/tranqui/rimeflows/blob/main/data/kuwabara-trajectories.png)

Run this to show the square root efficiency scaling for the Kuwabara flow field:

    python show-flow.py 0.15 -cm -ce data/efficiency_kuwabara_alpha\=0.15.csv

Which should generate this figure:

![Kuwabara particle trajectories](https://github.com/tranqui/rimeflows/blob/main/data/kuwabara-efficiency.png)

The various options available to customise output are shown in show-flow.py help page (i.e. by executing `python show-flow.py -h`).

The above commands make use of precalculated efficiency data which were generated for volume fractions alpha=0.1 and alpha=0.15 via the following commands:

    python efficiency.py 0.1 > data/efficiency_kuwabara_alpha=0.1.csv
    python efficiency.py 0.15 > data/efficiency_kuwabara_alpha=0.15.csv

These can be modified as needed for different values of alpha.

The flow fields implemented are:
* Toy stokes field: limiting flow field approaching a generic no-slip/incompressible surface
* Toy inviscid field: limiting flow field approaching a generic surface with slip boundary conditions
* Kuwabara
* Inviscid
* Hiemenz
* ?

### Bonus figures

Bonus figure, showing transition from exp scaling to normal (square root) scaling of capture efficiency as initial conditions are rescaled in the inviscid flow case:

    python inviscid-efficiency-vs-ics.py data/inviscid_eff_*

### Prepackaged data

Some of the figures above rely on precalculated data. This data was generated via:

    python efficiency.py 0.1 > data/efficiency_kuwabara_alpha=0.1.csv
    python efficiency.py 0.15 > data/efficiency_kuwabara_alpha=0.15.csv
    python efficiency.py -f stokes > data/efficiency_stokes.csv
    python efficiency.py -f shm > data/efficiency_shm.csv

### Out of date, need to update:

    for gamma in 1.2 1.8 2.0 2.5 3.0 3.5 5.0 10.0; do python efficiency.py --flow inviscid --scalefac $gamma --maxdeltastokes 1.0 --dump > data/efficiency_inviscid_gamma=$gamma.csv; done

Generating animated GIF showing change in streamlines as we move from limiting Stokes to inviscid case (SHM). First, create the streamlines for each flow field via:

    for eps in $(seq 0 0.01 1); do python on-axis-phase-portrait.py $eps; done

Or speed this up by doing it in parallel using xargs (replace $ncores with the number of CPU cores to use):

    seq 0 0.01 1 | xargs -P$ncores -I% python on-axis-phase-portrait.py %
Create the GIF:

    convert -delay 5 -loop 0 data/hybrid_streamlines_eps*.png -compress JPEG -quality 0 data/hybrid_streamlines_animated.gif


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

These programs are copyright &copy; 2022 Joshua Robinson and copyright
&copy; 2022 Patrick B Warren, as indicated in the individual files.
