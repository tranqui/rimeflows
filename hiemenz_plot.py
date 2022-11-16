#!/usr/bin/env python3

# hiemenz_plot.py -- plot or generate data for Hiemenz function.

# Copyright (c) 2022 Patrick B Warren <patrick.warren@stfc.ac.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Plot the Hiemenz function and its derivatives, or write out as a table

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utilities import Hiemenz

parser = argparse.ArgumentParser(description='plot Hiemenz function')
parser.add_argument('--hiemenz-cutoff', default=10.0, type=float, help='Hiemenz cutoff, default 10.0')
parser.add_argument('--hiemenz-bracket', default='1.22,1.24', help='Hiemenz bracket, default 1.22,1.24')
parser.add_argument('--hiemenz-tol', default=1.0e-12, type=float, help='Hiemenz tolerance, default 1.0e-12')
parser.add_argument('--linspace', default='0,20,80', help='linspace range for eta, default 0,20,80')
parser.add_argument('-o', '--output', help='write results as tab-separated values to a file')
args = parser.parse_args()

hiemenz = Hiemenz(cutoff=args.hiemenz_cutoff, init_bracket=eval('[' + args.hiemenz_bracket + ']'), tol=args.hiemenz_tol)

eta = np.linspace(*eval('[' + args.linspace + ']'))

phi, phip, phipp = hiemenz.func(eta)

df = pd.DataFrame({"eta":eta, "phi":phi, "phi'":phip, "phi''":phipp})

if args.output:

    df.set_index('eta', inplace=True) # avoid duplicating existing column

    with open(args.output, 'w') as f:
        print('#' + df.index.name + '\t' + '\t'.join(df.columns), file=f)
        print(df.to_csv(sep='\t', header=False), end='', file=f)

    print(f'Data written to {args.output}')

else: # plot phi / eta, phi', phi" as a function of eta

    df.index = df['eta'] # duplicates existing column but needed below
    df['phi / eta'] = df['phi'] / df['eta']
    df.at[0, 'phi/eta'] = 0 # because we know phi, phi' → 0 as eta → 0.
    df[["phi / eta", "phi'", "phi''"]].plot()
    plt.show()

