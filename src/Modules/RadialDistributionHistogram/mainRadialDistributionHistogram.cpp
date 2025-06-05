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
// Created by galavizp on 5/06/2025.
//

#include "mainRadialDistributionHistogram.h"
#include "trajectoryReader.h"
#include "logger.h"
#include <boost/histogram.hpp>

namespace mdtools {

    void mainRadialDistributionHistogram(const radial_distribution_histogram_options_t &radial_distribution_histogram,
                                         const io_options_t &io_options, simulation_options_t simulation_options) {

        auto trajectory = trajectoryReader(io_options.trajectory_input_file, io_options.coordinates_input_file).get(
                simulation_options.time_step, simulation_options.start_iteration, simulation_options.delta_iteration,
                simulation_options.end_iteration);

        if (trajectory.empty()) {
            LOGGER.error << "RadialDistributionHistogram failed" << std::endl;
            return;
        }

        auto number_of_frames = trajectory[0].position_x.size();
        LOGGER.info << "Reading done. Number of frames: " << number_of_frames << std::endl;

        std::map<int, boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>> histograms{};

        auto center = str2center[radial_distribution_histogram.center];
        double total_mass = 0.0;
        for (auto &atom: trajectory) {
            if (histograms.find(atom.atom_type) == histograms.end()) {
                histograms[atom.atom_type] = boost::histogram::make_histogram(
                        boost::histogram::axis::regular<>(radial_distribution_histogram.size,
                                                          radial_distribution_histogram.start,
                                                          radial_distribution_histogram.stop, "r"));
            }
            if (simulation_options.mass_map.find(atom.atom_type) == simulation_options.mass_map.end()) {
                LOGGER.error << "Unknown atom type:" << atom.atom_type << ". Define mass in simulation.atom_mass"
                             << std::endl;
                std::throw_with_nested(std::runtime_error("Fatal error"));
            }

            total_mass += simulation_options.mass_map[atom.atom_type];
        }

        double cx, cy, cz;
        for (int iteration = 0; iteration < number_of_frames; iteration++) {

            if (center == center_t::ORIGIN) {
                cx = 0.5 * (trajectory[0].lattice_origin_x[iteration] + trajectory[0].lattice_a[iteration]);
                cy = 0.5 * (trajectory[0].lattice_origin_y[iteration] + trajectory[0].lattice_b[iteration]);
                cz = 0.5 * (trajectory[0].lattice_origin_z[iteration] + trajectory[0].lattice_c[iteration]);
            } else {
                cx = 0;
                cy = 0;
                cz = 0;
                for (auto &atom: trajectory) {
                    auto mass = simulation_options.mass_map[atom.atom_type];
                    cx += mass * atom.position_x[iteration];
                    cy += mass * atom.position_y[iteration];
                    cz += mass * atom.position_z[iteration];
                }
                cx /= total_mass;
                cy /= total_mass;
                cz /= total_mass;
            }

            for (auto &atom: trajectory) {
                auto x = atom.position_x[iteration] - cx;
                auto y = atom.position_y[iteration] - cy;
                auto z = atom.position_z[iteration] - cz;
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