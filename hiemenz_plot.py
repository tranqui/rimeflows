#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
parser.add_argument('--linspace', default='0,20,80', help='linspace range for η, default 0,20,80')
parser.add_argument('-o', '--output', help='write results as tab-separated values to a file')
args = parser.parse_args()

hiemenz = Hiemenz(cutoff=args.hiemenz_cutoff, init_bracket=eval('[' + args.hiemenz_bracket + ']'), tol=args.hiemenz_tol)

η = np.linspace(*eval('[' + args.linspace + ']'))

φ, φp, φpp = hiemenz.func(η)

df = pd.DataFrame({"η":η, "φ":φ, "φ'":φp, "φ''":φpp})

if args.output:

    df.set_index('η', inplace=True) # avoid duplicating existing column

    with open(args.output, 'w') as f:
        print('#' + df.index.name + '\t' + '\t'.join(df.columns), file=f)
        print(df.to_csv(sep='\t', header=False), end='', file=f)

    print(f'Data written to {args.output}')

else: # plot φ / η, φ', φ" as a function of η

    df.index = df['η'] # duplicates existing column but needed below
    df['φ / η'] = df['φ'] / df['η']
    df.at[0, 'φ/η'] = 0 # because we know φ, φ' → 0 as η → 0.
    df[["φ / η", "φ'", "φ''"]].plot()
    plt.show()

