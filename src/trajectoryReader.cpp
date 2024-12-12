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
// Created by Pablo Galaviz on 11/4/2024.
//

#include "trajectoryReader.h"
#include <iostream>


namespace mdtools {
        trajectoryReader::trajectoryReader(const std::string& file_name) {


        auto format_flag = file_format(file_name);
        auto format = static_cast<format_t>(format_flag & format_t::FILE_TYPE);

        file.open(file_name, std::ios_base::in | std::ios_base::binary);
        if (format_flag & GZIP) { input_buffer.push(boost::iostreams::gzip_decompressor()); }
        input_buffer.push(file);
        //Convert stream buffer to istream
        input_stream = new std::istream(&input_buffer);
        input_buffer.set_auto_close(false);

        switch (format) {

            case format_t::LAMMPS:
            {getTrajectory = &trajectoryReader::getLammpsTrajectory;}
                break;
            case format_t::XDATCAR:
            {getTrajectory = &trajectoryReader::getXDATCARTrajectory;}
                break;

            case format_t::UNKNOWN:
            default:
                LOGGER.error << "Unknown input file format " << std::endl;
                exit(EINVAL);
        }

    }

    std::vector<atom_t> trajectoryReader::getLammpsTrajectory(double time_step) {

        std::string line;

        int number_of_atoms = 0;
        size_t box_coordinate = 0;
        state_t current_state = READING_NONE;
        std::vector<frame> trajectory;
        frame new_frame;
        while (std::getline(*input_stream, line)) {

            std::size_t item_found = line.find("ITEM");
            if (item_found != std::string::npos) {
                std::vector<std::string> fields;
                boost::algorithm::split(fields, line, boost::is_any_of(":"));
                if (fields[1].find("TIMESTEP") != std::string::npos) {
                    current_state = READING_STEP;
                    if (new_frame.time_step_id >= 0) {
                        trajectory.push_back(new_frame);
                        new_frame.reset();
                    }
                    continue;
                }
                if (fields[1].find("NUMBER OF ATOMS") != std::string::npos) {
                    current_state = READING_N_ATOMS;
                    continue;
                }
                if (fields[1].find("BOX BOUNDS") != std::string::npos) {
                    current_state = READING_BOX;
                    continue;
                }
                if (fields[1].find("ATOMS") != std::string::npos) {
                    current_state = READING_POSITION;
                    continue;
                }
                LOGGER.info << fields[1] << std::endl;
            } else {
                if (current_state == READING_STEP) {
                    new_frame.time_step_id = stoi(line);
                    current_state = READING_NONE;
                    continue;
                }
                if (current_state == READING_N_ATOMS) {
                    size_t noa = stoi(line);
                    if (number_of_atoms == 0) {
                        number_of_atoms = noa;
                    } else {
                        if (number_of_atoms != noa) {
                            LOGGER.warning << "Inconsistent number of atoms at frame:" << new_frame.time_step_id
                                           << std::endl;
                            LOGGER.warning << "Found:" << noa << " expecting:" << number_of_atoms << std::endl;
                            exit(-1);
                        }
                    }
                    new_frame.number_of_atoms = number_of_atoms;
                    new_frame.position_x.resize(number_of_atoms);
                    new_frame.position_y.resize(number_of_atoms);
                    new_frame.position_z.resize(number_of_atoms);
                    new_frame.atom_type.resize(number_of_atoms);
                    current_state = READING_NONE;
                    continue;
                }
                if (current_state == READING_BOX) {
                    char *pEnd;
                    new_frame.lattice[box_coordinate].minimum = strtod(line.c_str(), &pEnd);
                    new_frame.lattice[box_coordinate].maximum = strtod(pEnd, NULL);
                    box_coordinate = (box_coordinate + 1) % 3;
                    continue;
                }

                if (current_state == READING_POSITION) {
                    char *pEnd;
                    auto id = strtol(line.c_str(), &pEnd, 10) - 1;
                    if (id < 0 || id >= new_frame.number_of_atoms) {
                        LOGGER.warning << line << std::endl;
                        continue;
                    }
                    new_frame.atom_type[id] = static_cast<int>(strtol(pEnd, &pEnd, 10));
                    new_frame.position_x[id] = strtod(pEnd, &pEnd);
                    new_frame.position_y[id] = strtod(pEnd, &pEnd);
                    new_frame.position_z[id] = strtod(pEnd, &pEnd);
                    continue;
                }

                LOGGER.debug << line << std::endl;
            }

        }

        std::vector<atom_t> atom_trajectory = std::vector<atom_t>(trajectory[0].number_of_atoms);
        for (int atom_id = 0; atom_id < atom_trajectory.size(); atom_id++) {
            atom_trajectory[atom_id].time_step = time_step;
            atom_trajectory[atom_id].atom_type = trajectory[0].atom_type[atom_id];
            atom_trajectory[atom_id].position_x.resize(trajectory.size());
            atom_trajectory[atom_id].position_y.resize(trajectory.size());
            atom_trajectory[atom_id].position_z.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_a.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_b.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_c.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_origin_x.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_origin_y.resize(trajectory.size());
            atom_trajectory[atom_id].lattice_origin_z.resize(trajectory.size());
            atom_trajectory[atom_id].time.resize(trajectory.size());
            for (int trajectory_id = 0; trajectory_id < trajectory.size(); trajectory_id++) {
                atom_trajectory[atom_id].position_x[trajectory_id] = trajectory[trajectory_id].position_x[atom_id];
                atom_trajectory[atom_id].position_y[trajectory_id] = trajectory[trajectory_id].position_y[atom_id];
                atom_trajectory[atom_id].position_z[trajectory_id] = trajectory[trajectory_id].position_z[atom_id];
                atom_trajectory[atom_id].time[trajectory_id] = trajectory[trajectory_id].time_step_id * time_step;
                atom_trajectory[atom_id].lattice_origin_x[trajectory_id] = trajectory[trajectory_id].lattice[X].minimum;
                atom_trajectory[atom_id].lattice_origin_y[trajectory_id] = trajectory[trajectory_id].lattice[Y].minimum;
                atom_trajectory[atom_id].lattice_origin_z[trajectory_id] = trajectory[trajectory_id].lattice[Z].minimum;
                atom_trajectory[atom_id].lattice_a[trajectory_id] = trajectory[trajectory_id].lattice[X].maximum-trajectory[trajectory_id].lattice[X].minimum;
                atom_trajectory[atom_id].lattice_b[trajectory_id] = trajectory[trajectory_id].lattice[Y].maximum-trajectory[trajectory_id].lattice[Y].minimum;
                atom_trajectory[atom_id].lattice_c[trajectory_id] = trajectory[trajectory_id].lattice[Z].maximum-trajectory[trajectory_id].lattice[Z].minimum;

            }
            atom_trajectory[atom_id].calculate_velocity();
        }

        return atom_trajectory;

    }

    std::vector<atom_t> trajectoryReader::getXDATCARTrajectory(double time_step) {
        return std::vector<atom_t>();
    }


    trajectoryReader::~trajectoryReader() = default;
} // mdtools