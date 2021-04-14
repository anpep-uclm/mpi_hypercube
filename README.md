<div align="center">
  <img width="200" src="http://i.4pcdn.org/pol/1564264361687.png">
  <h3>mpi_hypercube</h3>
  <blockquote>Toy program implementing a hypercube-interconnect network topology using MPI.</blockquote>
</div>

## Build Instructions
```shell
$ git clone https://github.com/anpep/mpi_hypercube
$ cd mpi_hypercube
$ make
$ ./tools/generate_input.py 16 > input.dat
$ mpirun --oversubscribe -n 17 mpi_hypercube 4 input.dat
```

## Open-source license
```
mpi_hypercube -- Implements a hypercube-interconnect network topology using OpenMPI
Copyright (c) 2021 Ángel Pérez <angel@ttm.sh>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 2.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
```