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

#ifndef MDTOOLS_TRAJECTORYREADER_H
#define MDTOOLS_TRAJECTORYREADER_H

#include <string>
#include <valarray>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include "logger.h"

namespace mdtools {

    // First 4 bits is file ID, the fifth bit is compression flag.
    enum format_t : int {
        UNKNOWN = 0b0000'0000, LAMMPS = 0b0000'0001, XDATCAR = 0b0000'0010, LAMMPS_GZ = 0b0001'0001,
        GZIP = 0b0001'0000 , FILE_TYPE = 0b0000'1111,
    };

    enum axis_t : int {X = 0, Y = 1, Z = 2};

    inline format_t file_format(std::string name) {

        std::vector<std::string> name_vec;
        boost::algorithm::split(name_vec, name, boost::is_any_of("."));
        size_t name_size = name_vec.size();
        if (name_size > 1) {
            std::string ext1 = name_vec[name_size - 1];
            bool compressed = false;
            if (ext1 == "gzip" || ext1 == "gz") {
                if (name_size == 2) { return format_t::UNKNOWN; }
                ext1 = name_vec[name_size - 2];
                compressed = true;
            }

            if (ext1 == "lammpstrj") {
                if (compressed) { return format_t::LAMMPS_GZ; }
                return format_t::LAMMPS;
            }
        }
        else{
            if(name_size ==1){
                if (name == "XDATCAR") {
                    return format_t::XDATCAR;
                }

            }
        }

        return format_t::UNKNOWN;
    }

    enum state_t : int {
        READING_NONE = 0, READING_STEP, READING_N_ATOMS, READING_BOX, READING_POSITION
    };


    struct box {
        double minimum = 0;
        double maximum = 1;
    };

    struct frame {
        int time_step_id = -1;
        size_t number_of_atoms = 0;
        std::vector<box> lattice = std::vector<box>(3);
        std::valarray<double> position_x;
        std::valarray<double> position_y;
        std::valarray<double> position_z;
        std::valarray<int> atom_type;

        void reset() {
            time_step_id = -1;
            number_of_atoms = 0;
            lattice = std::vector<box>(3);
            position_x.apply([](double x) -> double { return 0.0; });
            position_y.apply([](double x) -> double { return 0.0; });
            position_z.apply([](double x) -> double { return 0.0; });
            atom_type.apply([](int x) -> int { return 0; });
        }
    };


    static void periodic_boundary_correction(std::valarray<double> &position, const double &dt, std::valarray<double> &velocity){
        auto p0 = position[0];
        double idt = 1.0/dt;
//        LOGGER.debug << abs(velocity).max() << " -> ";

        for(int i=0; i < position.size()-1;i++){
            auto pi=position[i];
            auto pip1=position[i+1];
            if(abs(pi-p0) > 0.5){
                pi= pi>p0 ? pi-1 : pi+1;
                position[i]=pi;
            }
            if(abs(pip1-p0) > 0.5){
                pip1= pip1>p0 ? pip1-1 : pip1+1;
                position[i+1]=pip1;
            }

            velocity[i+1]=(pip1-pi)*idt;
        }
        velocity[0]=velocity[1];
//        LOGGER.debug << abs(velocity).max() << std::endl;
    }

    struct atom_t {
        double time_step=0;
        int atom_type=0;
        std::valarray<double> position_x;
        std::valarray<double> position_y;
        std::valarray<double> position_z;

        std::valarray<double> velocity_x;
        std::valarray<double> velocity_y;
        std::valarray<double> velocity_z;

        std::valarray<double> time;
        std::valarray<double> lattice_a;
        std::valarray<double> lattice_b;
        std::valarray<double> lattice_c;

        std::valarray<double> lattice_origin_x;
        std::valarray<double> lattice_origin_y;
        std::valarray<double> lattice_origin_z;


        std::string serialize(){
            std::stringstream  result;
            result << "{";
            result << "\"time step\" : " << time_step << ",";
            result << "\"atom_t type\" : " << atom_type << ",";
            result << "\"position x\" : [ ";
            for(auto item : position_x){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"position y\" : [ ";
            for(auto item : position_y){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"position z\" : [ ";
            for(auto item : position_z){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"velocity x\" : [ ";
            for(auto item : velocity_x){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"velocity y\" : [ ";
            for(auto item : velocity_y){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"velocity z\" : [ ";
            for(auto item : velocity_z){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "],";
            result << "\"time\" : [ ";
            for(auto item : time){result << item << ",";}
            result.seekp(-1, std::ios_base::end);
            result << "]";
            result << "}";
            return result.str();
        }

        void calculate_velocity(){
            double idt= 1.0/time_step;
            velocity_x = idt*(position_x-position_x.shift(-1));
            velocity_x[0]=velocity_x[1];
            velocity_y = idt*(position_y-position_y.shift(-1));
            velocity_y[0]=velocity_y[1];
            velocity_z = idt*(position_z-position_z.shift(-1));
            velocity_z[0]=velocity_z[1];

            auto max_x=abs(velocity_x).max()*time_step;
            if(max_x > 0.5){periodic_boundary_correction(position_x,time_step,velocity_x);}

            auto max_y=abs(velocity_y).max()*time_step;
            if(max_y > 0.5){periodic_boundary_correction(position_y,time_step,velocity_y);}

            auto max_z=abs(velocity_z).max()*time_step;
            if(max_z > 0.5){periodic_boundary_correction(position_z,time_step,velocity_z);}

            position_x = lattice_origin_x+lattice_a*position_x;
            position_y = lattice_origin_y+lattice_b*position_y;
            position_z = lattice_origin_z+lattice_c*position_z;

            velocity_x *= lattice_a;
            velocity_y *= lattice_b;
            velocity_z *= lattice_c;

        }

    };

    class trajectoryReader {

        std::istream *input_stream;
        std::ifstream file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> input_buffer;


        std::vector<atom_t> (trajectoryReader::*getTrajectory)(double time_step);
        std::vector<atom_t> getLammpsTrajectory(double time_step);
        std::vector<atom_t> getXDATCARTrajectory(double time_step);

    public:
        explicit trajectoryReader(const std::string& file_name);

        virtual ~trajectoryReader();

        inline std::vector<atom_t> get(double time_step) {return (this->*getTrajectory)(time_step);}

    };

} // mdtools

#endif //MDTOOLS_TRAJECTORYREADER_H
