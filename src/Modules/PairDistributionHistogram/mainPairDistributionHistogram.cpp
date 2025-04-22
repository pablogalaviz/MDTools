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

#include "mainPairDistributionHistogram.h"
#include "trajectoryReader.h"
#include "logger.h"
#include <boost/histogram.hpp>

namespace mdtools {

    void mainPairDistributionHistogram(const pair_distribution_histogram_options_t& pair_distribution_histogram,const io_options_t& io_options, simulation_options_t simulation_options){

        auto trajectory = trajectoryReader(io_options.trajectory_input_file,io_options.coordinates_input_file).get(simulation_options.time_step, simulation_options.start_iteration, simulation_options.end_iteration);

        if(trajectory.empty()){
            LOGGER.error << "PairDistributionHistogram failed"<< std::endl;
            return;
        }

        auto number_of_frames = trajectory[0].position_x.size();
        LOGGER.info << "Reading done. Number of frames: " << number_of_frames << std::endl;
        auto number_of_atoms = trajectory.size();

        std::map<int, boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>> histograms{};

#pragma omp parallel for schedule(static)
        for (int i = 0; i < number_of_atoms - 1; ++i) {
            auto &atom_i = trajectory[i];
            auto histogram_iter = histograms.find(atom_i.atom_type);

            // Construct a new histogram when atom_i.atom_type is not already present
#pragma omp critical
            {
                if (histogram_iter == histograms.end()) {
                    histogram_iter = histograms.insert({atom_i.atom_type, boost::histogram::make_histogram(
                            boost::histogram::axis::regular<>(pair_distribution_histogram.size,
                                                              pair_distribution_histogram.start,
                                                              pair_distribution_histogram.stop, "r"))}).first;
                }
            }

#pragma omp parallel for schedule(static)
            for (int j = i + 1; j < number_of_atoms; ++j) {
                auto &atom_j = trajectory[j];
                if (atom_i.atom_type != atom_j.atom_type) { continue; }

                double x = atom_i.mean_position_x - atom_j.mean_position_x;
                double y = atom_i.mean_position_y - atom_j.mean_position_y;
                double z = atom_i.mean_position_z - atom_j.mean_position_z;
                auto r = sqrt(x * x + y * y + z * z);

                // Use iterator to access histogram
#pragma omp critical
                {
                    histogram_iter->second(r);
                }
            }
        }

        for(auto item : histograms) {
            std::ofstream file;
            file.open(io_options.output_path + "/hist_" + std::to_string(item.first) + ".csv", std::ios_base::out);
            for (auto &&x: indexed(item.second, boost::histogram::coverage::all)) {
                file << 0.5 * (x.bin().lower() + x.bin().upper()) << "," << *x << std::endl;
            }
            file.close();
        }
    }

} // mdtools