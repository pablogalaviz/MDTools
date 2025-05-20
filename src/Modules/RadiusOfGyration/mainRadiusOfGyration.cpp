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

#include "mainRadiusOfGyration.h"
#include "trajectoryReader.h"
#include "logger.h"
#include <gsl/gsl_math.h>

namespace mdtools {

    void mainRadiusOfGyration(radius_of_gyration_options_t radius_of_gyration_options,
                              const io_options_t &io_options,
                              simulation_options_t simulation_options) {

        auto trajectory = trajectoryReader(io_options.trajectory_input_file, io_options.coordinates_input_file).get(
                simulation_options.time_step, simulation_options.start_iteration, simulation_options.delta_iteration,
                simulation_options.end_iteration);

        if (trajectory.empty()) {
            LOGGER.error << "RadiusOfGyration failed" << std::endl;
            return;
        }

        auto number_of_frames = trajectory[0].position_x.size();
        LOGGER.info << "Reading done. Number of frames: " << number_of_frames << std::endl;
        auto number_of_atoms = trajectory.size();


        double total_mass = 0.0;
        for (auto &atom: trajectory) {
            if (simulation_options.mass_map.find(atom.atom_type) == simulation_options.mass_map.end()) {
                LOGGER.error << "Unknown atom type:" << atom.atom_type << ". Define mass in simulation.atom_mass"
                             << std::endl;
                std::throw_with_nested(std::runtime_error("Fatal error"));
            }

            total_mass += simulation_options.mass_map[atom.atom_type];
        }

        std::ofstream file;
        file.open(io_options.output_path + "/rog.csv", std::ios_base::out);
        file << "Time (ps),Radius (nm)" << std::endl;
        for (int iteration = 0; iteration < number_of_frames; iteration++) {
            double center_of_mass_x = 0;
            double center_of_mass_y = 0;
            double center_of_mass_z = 0;
            for (auto &atom: trajectory) {
                auto mass = simulation_options.mass_map[atom.atom_type];
                center_of_mass_x += mass * atom.position_x[iteration];
                center_of_mass_y += mass * atom.position_y[iteration];
                center_of_mass_z += mass * atom.position_z[iteration];
            }

            center_of_mass_x /= total_mass;
            center_of_mass_y /= total_mass;
            center_of_mass_z /= total_mass;

            double radius_of_gyrate = 0;

            double max_x = 0;
            for (auto &atom: trajectory) {
                auto mass = simulation_options.mass_map[atom.atom_type];
                radius_of_gyrate += mass * (gsl_pow_2(atom.position_x[iteration] - center_of_mass_x)
                                            + gsl_pow_2(atom.position_y[iteration] - center_of_mass_y)
                                            + gsl_pow_2(atom.position_z[iteration] - center_of_mass_z));
                max_x = std::max(max_x, atom.position_x[iteration]);
            }
            radius_of_gyrate /= total_mass;

            file << trajectory[0].time[iteration] << "," << std::sqrt(radius_of_gyrate) << std::endl;
//            LOGGER.debug << trajectory[0].time[iteration] << "," << std::sqrt(radius_of_gyrate) << std::endl;

        }


    }

}
