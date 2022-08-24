#!/usr/bin/bash

# hiemenz.sh -- automate generation of data files

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

# The code generates commands which can be executed by piping the
# output into /usr/bin/bash.  If necessary this can be filtered, eg:

# ./hiemenz.sh | /usr/bin/bash
# ./hiemenz.sh | grep delta | /usr/bin/bash
# ./hiemenz.sh | grep reynolds | /usr/bin/bash

# Boundary layer thickness control (original case)

stokes=0.2,1.0,0.01
for x in 000 001 002 005 010 020
do
    echo "./hiemenz.py --delta=$x/100 --stokes=$stokes -o hiemenz$x.dat"
done

stokes=0.2,1.5,0.01
x=050
echo "./hiemenz.py --delta=$x/100 --stokes=$stokes -o hiemenz$x.dat"

stokes=0.5,2.0,0.01
x=100
echo "./hiemenz.py --delta=$x/100 --stokes=$stokes -o hiemenz$x.dat"

# Reynolds number control

stokes=0.25,1.0,0.01
x=inft
echo "./hiemenz.py --reynolds=$x --stokes=$stokes -o reynolds$x.dat"

stokes=0.20,1.0,0.01
for x in 0005 0010 0050 0100 0200 0500 1000 2000 5000
do
    echo "./hiemenz.py --reynolds=$x --stokes=$stokes -o reynolds$x.dat"
done

stokes=0.20,1.5,0.01
x=0002
echo "./hiemenz.py --reynolds=$x --stokes=$stokes -o reynolds$x.dat"

stokes=0.50,2.0,0.01
x=0001
echo "./hiemenz.py --reynolds=$x --stokes=$stokes -o reynolds$x.dat"
