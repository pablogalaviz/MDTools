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

#include "mainAxialDistributionHistogram.h"
#include "trajectoryReader.h"
#include "logger.h"
#include <boost/histogram.hpp>

namespace mdtools {

    void mainAxialDistributionHistogram(axial_distribution_histogram_options_t axial_distribution_histogram,const io_options_t& io_options, simulation_options_t simulation_options){

        auto trajectory = trajectoryReader(io_options.input_file).get(simulation_options.time_step);

        LOGGER.info << "reading done " << std::endl;

        auto n = trajectory[0].position_x.size();

        using namespace boost::histogram;
        auto histogram1 = make_histogram(axis::regular<>(100, 14, 27, "r"));

        for (auto &atom: trajectory) {
            if(atom.atom_type!=4){continue;}
            for (int i = 0; i < n; i++) {
                auto cy = 0.5*(atom.lattice_origin_y[i] + atom.lattice_b[i]);
                auto cz = 0.5*(atom.lattice_origin_z[i] + atom.lattice_c[i]);
                auto y = atom.position_y[i] - cy;
                auto z = atom.position_z[i] - cz;
                auto r = sqrt(y*y+z*z);
                histogram1(r);
            }
        }

        std::ofstream file;
        file.open(io_options.output_path+"/hist_H2O.csv", std::ios_base::out );

        for (auto&& x : indexed(histogram1, coverage::all)) {
            file << 0.5*(x.bin().lower() + x.bin().upper()) << "," << *x << std::endl;
        }

    }

} // mdtools