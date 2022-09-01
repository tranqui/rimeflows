## Calculations of point particle capture efficiency in model 2d flows

The main output currently is the two figures for our upcoming manuscript.
Run via:

    python fig1-schematics.py

and

    python fig2-phase-portraits.py

And these should output the figures as PDFs and PNGs in the working directory.

Bonus figure, showing transition from exp scaling to normal (square root) scaling of capture efficiency as initial conditions are rescaled in the inviscid flow case:

    python inviscid-efficiency-vs-ics.py data/inviscid_eff_*

The prepackaged data was previously generated via:

    python efficiency.py 0.1 > data/efficiency_kuwabara_alpha=0.1.csv
    python efficiency.py -f stokes > data/efficiency_stokes.csv
    python efficiency.py -f shm > data/efficiency_shm.csv

Generate 3d surface for separatrix in our model toy problem (requires povray):

    python separatrix3d.py --povray > data/separatrix3d.inc
    python fig3-backdrop.py
    (cd data && povray +Q11 +A +AM2 +R5 +ua scene.pov)

And the resulting scene.png file should be generated inside the data folder.

Out of date, need to update:

    for gamma in 1.2 1.8 2.0 2.5 3.0 3.5 5.0 10.0; do python efficiency.py --flow inviscid --scalefac $gamma --maxdeltastokes 1.0 --dump > data/efficiency_inviscid_gamma=$gamma.csv; done

Generating animated GIF showing change in streamlines as we move from limiting Stokes to inviscid case (SHM). First, create the streamlines for each flow field via:

for eps in $(seq 0 0.002 0.12); do python fig3-backdrop.py $eps; done

Or speed this up by doing it in parallel using xargs (replace $ncores with the number of CPU cores to use):

    seq 0 0.002 0.12 | xargs -P$ncores -I% python fig3-backdrop.py %
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
