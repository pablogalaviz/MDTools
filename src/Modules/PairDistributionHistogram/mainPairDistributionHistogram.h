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
// Created by galavizp on 12/12/2024.
//

#ifndef MDTOOLS_MAINPAIRDISTRIBUTIONHISTOGRAM_H
#define MDTOOLS_MAINPAIRDISTRIBUTIONHISTOGRAM_H

#include "parameters.h"

namespace mdtools {

    void mainPairDistributionHistogram(const pair_distribution_histogram_options_t& pair_distribution_histogram,const io_options_t& io_options, simulation_options_t simulation_options);

} // mdtools

#endif //MDTOOLS_MAINPAIRDISTRIBUTIONHISTOGRAM_H
