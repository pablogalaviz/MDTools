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

#include "io.h"

namespace mdtools {

    /**
 * Shows the current arguments and their corresponding values.
 * @param vm is a boost variable_map object.
 */
    void show_options(boost::program_options::variables_map vm) {
        LOGGER.info << "Parameters :" << std::endl;
        LOGGER.info << "..............................................................." << std::endl;
        for (auto &item: vm) {
            std::stringstream ss;
            ss << std::setw(45) << std::left << item.first << " : ";
            boost::any value = item.second.value();
            if (is_type<std::string>(value)) { ss << item.second.as<std::string>(); }
            else {
                if (is_type<int>(value)) { ss << item.second.as<int>(); }
                else {
                    if (is_type<unsigned int>(value)) { ss << item.second.as<unsigned int>(); }
                    else {
                        if (is_type<double>(value)) { ss << item.second.as<double>(); }
                        else {
                            if (is_type<bool>(value)) { ss << (item.second.as<bool>() ? "true" : "false"); }
                            else {
                                if (is_type<size_t>(value)) { ss << item.second.as<size_t>(); }
                            }
                        }
                    }
                }
            }
            LOGGER.info << ss.str() << std::endl;
        }
        LOGGER.info << "..............................................................." << std::endl;
    }

/**
 * Show and initialize message
 * @param code_name
 */
    void initialize(const std::string &code_name) {
        LOGGER.info << "-----------------------------------------------------------------" << std::endl;
        LOGGER.info << R"(     ___      .__   __.      _______.___________.  ______  )" << std::endl;
        LOGGER.info << R"(    /   \     |  \ |  |     /       |           | /  __  \ )" << std::endl;
        LOGGER.info << R"(   /  ^  \    |   \|  |    |   (----`---|  |----`|  |  |  |)" << std::endl;
        LOGGER.info << R"(  /  /_\  \   |  . `  |     \   \       |  |     |  |  |  |)" << std::endl;
        LOGGER.info << R"( /  _____  \  |  |\   | .----)   |      |  |     |  `--'  |)" << std::endl;
        LOGGER.info << R"(/__/     \__\ |__| \__| |_______/       |__|      \______/ )" << std::endl;
        LOGGER.info << "" << std::endl;
        LOGGER.info << "---- Australian Nuclear Science and Technology Organisation -----" << std::endl;
        LOGGER.info << "Nuclear science and technology for the benefit of all Australians" << std::endl;
        LOGGER.info << "" << std::endl;
        LOGGER.info << " ===============================" << std::endl;
        LOGGER.info << " " << code_name << std::endl;
        LOGGER.info << " Author: Pablo Galaviz             " << std::endl;
        LOGGER.info << " galavizp@ansto.gov.au              " << std::endl;
        LOGGER.info << " Build: " << __DATE__ << " " << __TIME__ << std::endl;
        LOGGER.info << " ===============================" << std::endl;
        LOGGER.info << "" << std::endl;
    }

/**
 * Returns a string in time format
 * @param value time in specified units
 * @param unit the unit of the time
 * @return
 */
    std::string get_time_str(long value, const std::string &unit) {

        std::stringstream result;

        if (value > 0) {
            result << value << " " << unit;
            if (value > 1)
                result << "s ";
            else
                result << " ";
        }
        return result.str();
    }

/**
 * Shows a good bye message and performance metrics
 * @param start
 */
    void finalize(std::chrono::system_clock::time_point start) {
        auto stop = std::chrono::system_clock::now();

        auto delta_hours = std::chrono::duration_cast<std::chrono::hours>(stop - start).count();
        auto delta_minutes = std::chrono::duration_cast<std::chrono::minutes>(stop - start).count();
        auto delta_seconds = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
        auto delta_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

        int days = int(delta_hours / 24);
        auto hours = delta_hours - days * 24;
        auto minutes = delta_minutes - delta_hours * 60;
        auto seconds = delta_seconds - delta_minutes * 60;
        auto milliseconds = delta_milliseconds - delta_seconds * 1000;

        LOGGER.info << "Finished in "
                    << get_time_str(days, "day")
                    << get_time_str(hours, "hour")
                    << get_time_str(minutes, "minute")
                    << get_time_str(seconds, "second")
                    << milliseconds << " milliseconds " << std::endl;

        LOGGER.info << "All done! " << std::endl;

    }



}



