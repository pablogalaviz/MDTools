// MDTools
//     Copyright (C)  2024  Pablo Galaviz
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <https://www.gnu.org/licenses/>.

//
// Created by Pablo Galaviz on 8/4/2024.
//

#ifndef MDTOOLS_MAINPHONONDOS_H
#define MDTOOLS_MAINPHONONDOS_H

#include "parameters.h"

namespace mdtools {
    void mainPhononDOS(phonon_dos_options_t phonon_dos_options,
                       const io_options_t& io_options,
                       const simulation_options_t& simulation_options
                       );
}


#endif //MDTOOLS_MAINPHONONDOS_H
