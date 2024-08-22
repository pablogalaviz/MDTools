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

#ifndef MDTOOLS_LOGGER_H
#define MDTOOLS_LOGGER_H

#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#define LOGGER Logger::getInstance()

namespace mdtools {

    enum class log_t : int {
        DBG = 0, INFO = 1, WARNING = 2, ERROR = 3
    };

    class Logger {
    public :

        static Logger &getInstance() {
            static Logger instance;
            return instance;
        }

    private :

        Logger() {};

        class LoggerStream {

            bool m_new_line;
            log_t m_level;
            log_t m_current;
            std::string m_label;

            std::map<log_t, std::string> color_map = {
                    {log_t::DBG,     "\033[1;32m"},
                    {log_t::INFO,    "\033[1;94m"},
                    {log_t::WARNING, "\033[1;33m"},
                    {log_t::ERROR,   "\033[1;31m"}};

            template<class T>
            std::stringstream get_formatted_stream(const T out, bool console) {

                std::stringstream ss;
                if (m_new_line) {
                    std::time_t t = std::time(nullptr);
                    std::tm tm = *std::localtime(&t);
                    ss << (console ? color_map[m_current] : "")
                       << m_label
                       << std::put_time(&tm, "[%F %T] | ")
                       << (console ? "\033[0m" : "");
                }
                ss << out;

                return ss;
            }


        public:

            std::ofstream file;

            LoggerStream() = default;

            ~LoggerStream() = default;

            void init(const log_t level,
                      const std::string& log_file,
                      const log_t current,
                      const std::string& label) {
                m_new_line = true;
                m_current = current;
                m_label = label;
                m_level = level;
                file.open(log_file, std::ios::app);
            };

            void close() {
                if (file.is_open()) {
                    file.close();
                }
            }

            template<class T>
            LoggerStream &operator<<(const T out) {
                if (m_level <= m_current) {
                    file << get_formatted_stream(out, false).str();
                    file.flush();
                        if (m_current == log_t::ERROR) { std::cerr << get_formatted_stream(out, true).str(); }
                        else {
                            std::cout << get_formatted_stream(out, true).str();
                        }
                    if (m_new_line) { m_new_line = false; }
                }
                return *this;
            }

            LoggerStream &operator<<(std::ostream &(*pfun)(std::ostream &)) {
                if (m_level <= m_current) {
                    pfun(file);
                        if (m_current == log_t::ERROR) { pfun(std::cerr); }
                        else { pfun(std::cout); }
                    m_new_line = true;
                }
                return *this;
            }

        };

        int m_id{};

    public:

        Logger(Logger const &) = delete;

        void operator=(Logger const &) = delete;

        inline void init(const log_t level, const std::string &log_file) {
            debug.init(level, log_file, log_t::DBG, "Debug   ");
            info.init(level, log_file, log_t::INFO, "Info    ");
            warning.init(level, log_file, log_t::WARNING, "Warning ");
            error.init(level, log_file, log_t::ERROR, "Error   ");
        }

        LoggerStream debug;
        LoggerStream info;
        LoggerStream warning;
        LoggerStream error;

        void close() {
            error.close();
        }
    };


}

#endif //MDTOOLS_LOGGER_H
