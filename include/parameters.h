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
        PhononDOS = 1, DynamicStructureFactor
    };

    /// Map a string argument to a task enumerator.
    __attribute__((unused)) static std::map<std::string, task_t> str2task{
            {"PhononDOS",     task_t::PhononDOS},
            {"DynamicStructureFactor", task_t::DynamicStructureFactor}
    };


    struct io_options_t {

        bool backup;
        std::string output_path;
        std::string input_file;
        int progress = 0;

        void validate() const {
            open_file(input_file,"Input file is missing");
        }
    };

    struct simulation_options_t {

        double time_step=0;

        void validate() const {
            if (time_step <= 0) {
                std::throw_with_nested(
                        std::runtime_error("Negative or zero time_step"));
            }
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


}

#endif //MDTOOLS_PARAMETERS_H
