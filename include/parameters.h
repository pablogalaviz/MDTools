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

#ifndef MDTOOLS_PARAMETERS_H
#define MDTOOLS_PARAMETERS_H

#include <map>
#include <string>
#include "io.h"

namespace mdtools {
    /**
 * Unified header for each task options
 */

    /// Defines an enumerator for the tasks
    enum class task_t : int {
        PhononDOS = 1, DynamicStructureFactor, AxialDistributionHistogram
    };

    /// Map a string argument to a task enumerator.
    __attribute__((unused)) static std::map<std::string, task_t> str2task{
            {"PhononDOS",     task_t::PhononDOS},
            {"DynamicStructureFactor", task_t::DynamicStructureFactor},
            {"AxialDistributionHistogram", task_t::AxialDistributionHistogram}
    };

    /// Defines an enumerator for the axis
    enum class axis_t : int {
        X = 0, Y, Z
    };

    /// Map a string argument to a task enumerator.
    __attribute__((unused)) static std::map<std::string, axis_t> str2axis{
            {"X",     axis_t::X},
            {"x", axis_t::X},
            {"Y", axis_t::Y},
            {"y", axis_t::Y},
            {"Z", axis_t::Z},
            {"z", axis_t::Z},
    };

    struct io_options_t {

        bool backup;
        std::string output_path;
        std::string trajectory_input_file;
        std::string coordinates_input_file;
        int progress = 0;

        void validate() const {
            open_file(trajectory_input_file, "Trajectory input file is missing");
        }
    };

    struct simulation_options_t {

        double time_step=0;
        int start_iteration =0;
        int end_iteration = -1;

        void validate() {
            if (time_step <= 0) {
                std::throw_with_nested(
                        std::runtime_error("Negative or zero time_step"));
            }
            if (start_iteration < 0) {
                std::throw_with_nested(
                        std::runtime_error("Negative or zero start iteration"));
            }

            if (end_iteration <= start_iteration){end_iteration = -1;}

        }
    };


    struct phonon_dos_options_t {

        double sigma=0;

        void validate() const {

            if (sigma <= 0) {
                std::throw_with_nested(
                        std::runtime_error("Negative or zero sigma"));
            }

        }

    };

    struct dynamic_structure_factor_options_t {

        double val=0;

        void validate() const {

            if (val <= 0) {
                std::throw_with_nested(
                        std::runtime_error("Negative or zero val"));
            }


        }

    };

    struct axial_distribution_histogram_options_t {

        std::string axis="x";
        double start=0;
        double stop=1;
        int size=100;

        void validate() const {

            if (str2axis.find(axis) == str2axis.end()) {
                std::throw_with_nested(
                        std::runtime_error("axis should be 0,1 or 2"));
            }


        }

    };

}

#endif //MDTOOLS_PARAMETERS_H
