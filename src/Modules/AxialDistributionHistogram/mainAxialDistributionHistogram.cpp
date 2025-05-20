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

    void mainAxialDistributionHistogram(const axial_distribution_histogram_options_t &axial_distribution_histogram,
                                        const io_options_t &io_options, simulation_options_t simulation_options) {

        auto trajectory = trajectoryReader(io_options.trajectory_input_file, io_options.coordinates_input_file).get(
                simulation_options.time_step, simulation_options.start_iteration, simulation_options.delta_iteration,
                simulation_options.end_iteration);

        if (trajectory.empty()) {
            LOGGER.error << "AxialDistributionHistogram failed" << std::endl;
            return;
        }

        auto n = trajectory[0].position_x.size();
        LOGGER.info << "Reading done. Number of frames: " << n << std::endl;

        std::map<int, boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>> histograms{};

        auto axis = str2axis[axial_distribution_histogram.axis];
        for (auto &atom: trajectory) {
            if (histograms.find(atom.atom_type) == histograms.end()) {
                histograms[atom.atom_type] = boost::histogram::make_histogram(
                        boost::histogram::axis::regular<>(axial_distribution_histogram.size,
                                                          axial_distribution_histogram.start,
                                                          axial_distribution_histogram.stop, "r"));
            }
            for (int i = 0; i < n; i++) {
                auto cx = 0.5 * (atom.lattice_origin_x[i] + atom.lattice_a[i]);
                auto cy = 0.5 * (atom.lattice_origin_y[i] + atom.lattice_b[i]);
                auto cz = 0.5 * (atom.lattice_origin_z[i] + atom.lattice_c[i]);
                auto x = axis == axis_t::X ? 0 : atom.position_x[i] - cx;
                auto y = axis == axis_t::Y ? 0 : atom.position_y[i] - cy;
                auto z = axis == axis_t::Z ? 0 : atom.position_z[i] - cz;
                auto r = sqrt(x * x + y * y + z * z);
                histograms[atom.atom_type](r);
            }
        }

        for (auto item: histograms) {
            std::ofstream file;
            file.open(io_options.output_path + "/hist_" + std::to_string(item.first) + ".csv", std::ios_base::out);
            for (auto &&x: indexed(item.second, boost::histogram::coverage::all)) {
                file << 0.5 * (x.bin().lower() + x.bin().upper()) << "," << *x << std::endl;
            }
            file.close();
        }
    }

} // mdtools