// MDTools
//     Copyright (C)  2025  Pablo Galaviz
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
// Created by galavizp on 15/05/2025.
//

#ifndef MDTOOLS_MAINRADIUSOFGYRATION_H
#define MDTOOLS_MAINRADIUSOFGYRATION_H

#include "parameters.h"

namespace mdtools{
    void mainRadiusOfGyration(radius_of_gyration_options_t radius_of_gyration_options, const io_options_t& io_options, simulation_options_t simulation_options);
}

#endif //MDTOOLS_MAINRADIUSOFGYRATION_H
