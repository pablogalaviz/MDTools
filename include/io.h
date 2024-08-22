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

#ifndef MDTOOLS_IO_H
#define MDTOOLS_IO_H

#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include "logger.h"
#include <chrono>

namespace mdtools {
    inline void create_output_directory(std::string path, bool backup) {


            int sys_out;
            if (backup) {
                std::string rmdir = "rm -rf " + path + "_prev";
                sys_out = system(rmdir.c_str());
                if (sys_out) {
                    LOGGER.debug << "rm -rf return value " << sys_out << std::endl;
                }
                std::string mvdir = "mv -f " + path + " " + path + "_prev 2>/dev/null";
                sys_out = system(mvdir.c_str());
                if (sys_out) {
                    LOGGER.debug << "mv -f return value " << sys_out << std::endl;
                }
            } else {
                std::string rmdir = "rm -rf " + path;
                sys_out = system(rmdir.c_str());
            }

            std::string mkdir = "mkdir -p " + path;
            sys_out = system(mkdir.c_str());


    }

    inline void log_command(std::string path, const int ac, char *av[]) {


            std::string cmd = "echo '#!/bin/bash' > " + path + "/command.sh";
            int sys_out = system(cmd.c_str());

            cmd = "echo cd `pwd` >> " + path + "/command.sh";
            sys_out = system(cmd.c_str());

            std::stringstream param;

            for (int i = 0; i < ac; i++)

                param << av[i] << " ";

            cmd = "echo " + param.str() + " >> " + path + "/command.sh;";
            sys_out = system(cmd.c_str());


            cmd = "chmod +x " + path + "/command.sh";
            sys_out = system(cmd.c_str());


    }

    static std::ifstream open_file(const std::string &file_name, const std::string &error_message) {
        std::ifstream ifs(file_name.c_str());
        if (!ifs) {
            LOGGER.error << "No such file or directory: " << file_name << std::endl;
            LOGGER.error << error_message << std::endl;
            exit(ENOENT);
        }
        return ifs;
    }

    void show_options(boost::program_options::variables_map vm);

    void initialize(const std::string &code_name);

    std::string get_time_str(long value, const std::string& unit);

    void finalize(std::chrono::system_clock::time_point start);

    template<class T>
    bool is_type(const boost::any &operand) {
        return operand.type() == typeid(T);
    }




}

#endif //MDTOOLS_IO_H
