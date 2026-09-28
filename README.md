# MDTools

This repository contains the Molecular Dynamics Tools software developed by the Scientific Computing Team at the Australian Centre for Neutron Scattering (ACNS).

## Getting Started

git clone repository
```shell
git clone git@github.com:pablogalaviz/MDTools.git
cd MDTools
```


## Compilation & Installation

### Requirements

| Dependency | Minimum version | Purpose | Notes |
|---|---|---|---|
| CMake | 3.30 | Build system | `project(VERSION 3.30)` in `CMakeLists.txt` |
| C/C++ compiler | C++17 | Build | GCC >= 9, Clang >= 10, or Apple Clang >= 12 |
| Boost | 1.70 | `program_options`, `iostreams`, `filesystem`, `date_time` | Must be the **compiled** libraries, not headers-only |
| GSL | 2.x | Numerics | GNU Scientific Library |
| FFTW | 3.x | FFTs | Discovered via `pkg-config` (`fftw3`) |
| OpenMP | - | Parallelism | Ships with most compilers |
| pkg-config | - | FFTW discovery | Required, used by `pkg_check_modules` |

> **Note on Boost:** the namespaced targets `Boost::program_options`, `Boost::iostreams`,
> `Boost::filesystem`, and `Boost::date_time` are only created when the corresponding
> `COMPONENTS` are requested. You must install the **binary/compiled** Boost development
> packages (see below), not just the header-only package.

---

### Linux - Debian / Ubuntu

```bash  
sudo apt update  
sudo apt install -y \
    build-essential \
    cmake \
    pkg-config \
    libboost-program-options-dev \
    libboost-iostreams-dev \
    libboost-filesystem-dev \
    libboost-date-time-dev \
    libgsl-dev \
    libfftw3-dev \
    libomp-dev  

```


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