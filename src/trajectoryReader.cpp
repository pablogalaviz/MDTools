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
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include <iostream>

namespace mdtools {
        trajectoryReader::trajectoryReader(const std::string& fn,const std::string& coordinates_file_name) {

        file_name=fn;
        auto format_flag = file_format(file_name);
        auto format = static_cast<format_t>(format_flag & format_t::FILE_TYPE);

        switch (format) {

            case format_t::LAMMPS:
            {
                file.open(file_name, std::ios_base::in | std::ios_base::binary);
                if (format_flag & GZIP) { input_buffer.push(boost::iostreams::gzip_decompressor()); }
                input_buffer.push(file);
                //Convert stream buffer to istream
                input_stream = new std::istream(&input_buffer);
                input_buffer.set_auto_close(false);
                getTrajectory = &trajectoryReader::getLammpsTrajectory;
            }
                break;
            case format_t::XDATCAR:
            {getTrajectory = &trajectoryReader::getXDATCARTrajectory;}
                break;
            case format_t::XTC:
            {
                getTrajectory = &trajectoryReader::getXTCTrajectory;}
                break;
            case format_t::TRR:
            {
                auto format_flag = file_format(coordinates_file_name);
                auto format = static_cast<format_t>(format_flag & format_t::FILE_TYPE);
                if(format != format_t::GRO){
                    LOGGER.error << "Required gro coordinates file" << std::endl;
                    exit(EINVAL);
                }

                file.open(coordinates_file_name, std::ios_base::in);
                if (!file) {
                    LOGGER.error << "No such file or directory: " << coordinates_file_name << std::endl;
                    exit(ENOENT);
                }

                input_buffer.push(file);
                //Convert stream buffer to istream
                input_stream = new std::istream(&input_buffer);
                input_buffer.set_auto_close(false);
                getTrajectory = &trajectoryReader::getTRRTrajectory;}
                break;

            case format_t::UNKNOWN:
            default:
                LOGGER.error << "Unknown input file format " << std::endl;
                exit(EINVAL);
        }

    }

    std::vector<atom_t> trajectoryReader::getLammpsTrajectory(double time_step, int start_iteration,int delta_iteration ,int end_iteration) {

        std::string line;

        int number_of_atoms = 0;
        int frame_count = 0;
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
                    frame_count++;
                    if(frame_count <= start_iteration || frame_count%delta_iteration !=0 ){
                        current_state = SKIP_FRAME;
                    }
                    else
                    {
                        if(end_iteration > 0 && frame_count > end_iteration ){
                            break;
                        }
                        current_state = READING_STEP;
                        if (new_frame.time_step_id >= 0) {
                            trajectory.push_back(new_frame);
                            new_frame.reset();
                        }
                    }
                    continue;
                }
                if (current_state == SKIP_FRAME){
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
                if (current_state == SKIP_FRAME){
                    continue;
                }
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
                    new_frame.lattice[box_coordinate].minimum = strtod(line.c_str(), &pEnd)/10;
                    new_frame.lattice[box_coordinate].maximum = strtod(pEnd, NULL)/10;
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
            atom_trajectory[atom_id].calculate_means();
        }

        return atom_trajectory;

    }

    std::vector<atom_t> trajectoryReader::getXDATCARTrajectory(double time_step, int start_iteration, int delta_iteration, int end_iteration) {
        return std::vector<atom_t>();
    }

    std::vector<atom_t> trajectoryReader::getXTCTrajectory(double time_step, int start_iteration, int delta_iteration, int end_iteration) {
        return std::vector<atom_t>();
    }


    std::string trajectoryReader::removeNumbers(const std::string& input) {
        std::string result = input;
        result.erase(std::remove_if(result.begin(), result.end(), [](char c) { return !std::isalpha(c); }), result.end());
        return result;
    }

    std::vector<int> trajectoryReader::getAtomTypeFromGro(){
        std::vector<int> result;
        std::string line;
        int count=0;
        int atom_id=1;
        unsigned long number_of_atoms=0;
        std::map<std::string , int> map_atom_id;
        while (std::getline(*input_stream, line)) {
            if(count ==0){ count++; continue;}
            if(count ==1){ number_of_atoms=stoi(line) ;count++; continue;}
            std::vector<std::string> fields;
            boost::algorithm::split(fields, boost::trim_copy(line), boost::is_space(), boost::token_compress_on);
            if(fields.size() == 3){
                continue;
            }

            std::string atom_key=line.substr(9, 6);
            if(map_atom_id.find(atom_key) == map_atom_id.end()){
                map_atom_id[atom_key]=atom_id;
                atom_id++;
            }
            result.push_back(map_atom_id[atom_key]);
            count++;
        }

        if(result.size()!=number_of_atoms){
            LOGGER.warning << "Inconsistent number of atoms " << std::endl;
            LOGGER.warning << "Found:" << result.size() << " expecting:" << number_of_atoms << std::endl;
            exit(-1);
        }

        return result;
        }


    std::vector<atom_t> trajectoryReader::getTRRTrajectory(double time_step, int start_iteration, int delta_iteration, int end_iteration) {

         std::vector<int> atom_type = getAtomTypeFromGro();

        std::vector<atom_t> atom_trajectory{};
        int number_of_atoms;
        unsigned long number_of_frames;
        int64_t* offsets = nullptr;
        // Get number of atoms in the trajectory
        if (read_trr_header(file_name.c_str(), &number_of_atoms,&number_of_frames,&offsets) != exdrOK) {
            LOGGER.error << "Failed to read number of atoms from"<< file_name.c_str() << std::endl;
            exit(-1);
        }

        if (number_of_frames < start_iteration) {
            LOGGER.error << "Number of frames is smaller than simulation.start_iteration"<< std::endl;
            exit(-1);
        }

        if (end_iteration > 0 && number_of_frames < end_iteration) {
            LOGGER.warning << "Number of frames is smaller than simulation.end_iteration"<< std::endl;
        }

        // Allocate arrays for coordinates and velocities
        std::vector<float> coordinates(3 * number_of_atoms), velocity(3 * number_of_atoms);
        matrix box;
        int step;
        float time, lambda;
        XDRFILE* xdr = xdrfile_open(file_name.c_str(), "r");
        if (!xdr) {
            LOGGER.error << "Cannot open TRR file" << std::endl;
            return atom_trajectory;
        }

        number_of_frames=std::min(number_of_frames,static_cast<unsigned long>(end_iteration-start_iteration));

        atom_trajectory.resize(number_of_atoms);
        for (auto & atom : atom_trajectory) {
            atom.time_step = time_step;
            atom.position_x.resize(number_of_frames);
            atom.position_y.resize(number_of_frames);
            atom.position_z.resize(number_of_frames);
            atom.velocity_x.resize(number_of_frames);
            atom.velocity_y.resize(number_of_frames);
            atom.velocity_z.resize(number_of_frames);
            atom.lattice_a.resize(number_of_frames);
            atom.lattice_b.resize(number_of_frames);
            atom.lattice_c.resize(number_of_frames);
            atom.lattice_origin_x.resize(number_of_frames);
            atom.lattice_origin_y.resize(number_of_frames);
            atom.lattice_origin_z.resize(number_of_frames);
            atom.time.resize(number_of_frames);
        }
        unsigned long frame_id=0;
        int status;
        uint8_t flag = 0;
        while ((status = read_trr(xdr, number_of_atoms, &step, &time, &lambda,
                                  box,                         // simulation box (3x3 matrix)
                                  reinterpret_cast<rvec*>(coordinates.data()),  // coords array
                                  reinterpret_cast<rvec*>(velocity.data()),    // velocities array
                                  nullptr, &flag)) == exdrOK and frame_id < end_iteration) {
            if(frame_id >= start_iteration && frame_id%delta_iteration == 0){
                for(int atom_id=0; atom_id < number_of_atoms; atom_id++){
                    atom_trajectory[atom_id].position_x[frame_id-start_iteration]=coordinates[3 * atom_id];
                    atom_trajectory[atom_id].position_y[frame_id-start_iteration]=coordinates[3 * atom_id + 1];
                    atom_trajectory[atom_id].position_z[frame_id-start_iteration]=coordinates[3 * atom_id + 2];
                    atom_trajectory[atom_id].velocity_x[frame_id-start_iteration]=velocity[3 * atom_id];
                    atom_trajectory[atom_id].velocity_y[frame_id-start_iteration]=velocity[3 * atom_id + 1];
                    atom_trajectory[atom_id].velocity_z[frame_id-start_iteration]=velocity[3 * atom_id + 2];
                    atom_trajectory[atom_id].lattice_origin_x[frame_id-start_iteration]=0;
                    atom_trajectory[atom_id].lattice_origin_y[frame_id-start_iteration]=0;
                    atom_trajectory[atom_id].lattice_origin_z[frame_id-start_iteration]=0;
                    atom_trajectory[atom_id].lattice_a[frame_id-start_iteration]=box[0][0];
                    atom_trajectory[atom_id].lattice_b[frame_id-start_iteration]=box[1][1];
                    atom_trajectory[atom_id].lattice_c[frame_id-start_iteration]=box[2][2];
                    atom_trajectory[atom_id].time[frame_id-start_iteration]=step*time_step;
                    atom_trajectory[atom_id].atom_type=atom_type[atom_id];
                }
            }
            frame_id++;
        }

        for(auto &atom : atom_trajectory){
            atom.calculate_means();
        }

        xdrfile_close(xdr);

        return atom_trajectory;
    }

    trajectoryReader::~trajectoryReader() = default;
} // mdtools