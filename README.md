# MDTools

This repository contains the Molecular Dynamics Tools software developed by the Scientific Computing Team at the Australian Centre for Neutron Scattering (ACNS).

## Getting Started

git clone repository
```shell
git clone git@github.com:pablogalaviz/MDTools.git
cd MDTools
```
Install dependencies and compile using cmake tools
```shell
mkdir build
cd build
cmake ..
make 
make install
```

### Prerequisites

C++17 compiler, [Boost libraries](https://www.boost.org/), [GSL - GNU Scientific Library](https://www.gnu.org/software/gsl/) and [FFTW3](https://fftw.org/).

## History

First release July 2025

## Credits

Author: Pablo Galaviz

Contact: https://www.ansto.gov.au/scientific-computing


**Nuclear science and technology for the benefit of all Australians**

## License

MDTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

MDTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MDTools.  If not, see <http://www.gnu.org/licenses/>.