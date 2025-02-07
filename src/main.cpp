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

#include <iostream>
#include "io.h"
#include "parameters.h"
#include "Modules/DynamicStructureFactor/mainDynamicStructureFactor.h"
#include "Modules/PhononDOS/mainPhononDOS.h"
#include "Modules/AxialDistributionHistogram/mainAxialDistributionHistogram.h"

int main(const int ac, char *av[]) {

    int return_state = 0;
    try {

        auto start_time = std::chrono::system_clock::now();

        std::stringstream task_help;
        task_help << "Perform one of the following tasks: ";
        for (auto &item: mdtools::str2task) { task_help << item.first << " "; }

        std::string task{};
        std::string parameters;
        boost::program_options::options_description genericOptions(
                "Molecular Dynamics Analysis Tools.  \nOptions");
        genericOptions.add_options()
                ("debug,d", "Shows debug messages in log")
                ("help,h", "Shows a help message")
                ("task", boost::program_options::value<std::string>(&task), task_help.str().c_str())
                ("parameters,p", boost::program_options::value<std::string>(&parameters), "Parameters file")
                ("silent,s", "Shows only errors");

        std::string input_file;


        mdtools::io_options_t io_options;
        boost::program_options::options_description inputOutputOptions("Execution Options");
        inputOutputOptions.add_options()
                ("io.backup", boost::program_options::value<bool>(&io_options.backup)->default_value(true), "")
                ("io.output",
                 boost::program_options::value<std::string>(&io_options.output_path)->default_value("output"), "")
                ("io.input",
                 boost::program_options::value<std::string>(&io_options.input_file)->default_value("dump.lammpstrj"), "");


        mdtools::simulation_options_t simulation_options;
        boost::program_options::options_description simulationOptions("Simulation Options");
        simulationOptions.add_options()
                ("simulation.time_step",boost::program_options::value<double>(&simulation_options.time_step)->default_value(1), "Simulation time step in fs");


        mdtools::phonon_dos_options_t phonon_dos;
        boost::program_options::options_description phononDOSOptions("Phonon DOS Options");
        phononDOSOptions.add_options()
                ("phonon_dos.sigma",
                 boost::program_options::value<double>(&phonon_dos.sigma)->default_value(1), "");

        mdtools::dynamic_structure_factor_options_t dynamic_structure_factor;
        boost::program_options::options_description dynamicStructureFactorOptions("Dynamic Structure Factor Options");
        dynamicStructureFactorOptions.add_options()
                ("dynamic_structure_factor.val",
                 boost::program_options::value<double>(&dynamic_structure_factor.val)->default_value(1), "");

        mdtools::axial_distribution_histogram_options_t axial_distribution_histogram;
        boost::program_options::options_description axialDistributionHistogramOptions("Axial Distribution Histogram Options");
        axialDistributionHistogramOptions.add_options()
                ("axial_distribution_histogram.axis",
                 boost::program_options::value<std::string>(&axial_distribution_histogram.axis)->default_value("x"), "select the axis. Possible options [X,Y,Z,x,y,z]")
                ("axial_distribution_histogram.start",
                 boost::program_options::value<double>(&axial_distribution_histogram.start)->default_value(0), "Histogram start from the axis")
                ("axial_distribution_histogram.stop",
                 boost::program_options::value<double>(&axial_distribution_histogram.stop)->default_value(1), "Histogram stop from the axis")
                ("axial_distribution_histogram.size",
                 boost::program_options::value<int>(&axial_distribution_histogram.size)->default_value(100), "Histogram number of bins");



        boost::program_options::positional_options_description positional;
        positional.add("task", 1);

        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions)
                .add(inputOutputOptions)
                .add(simulationOptions)
                .add(phononDOSOptions)
                .add(dynamicStructureFactorOptions)
                .add(axialDistributionHistogramOptions)
                ;

        boost::program_options::options_description configFileOptions;
        configFileOptions.add(inputOutputOptions)
                .add(simulationOptions)
                .add(phononDOSOptions)
                .add(dynamicStructureFactorOptions)
                .add(axialDistributionHistogramOptions)
                ;

        boost::program_options::variables_map vm;
        boost::program_options::store(
                boost::program_options::command_line_parser(ac, av).options(cmdlineOptions).positional(
                        positional).run(), vm);
        boost::program_options::notify(vm);

        if (vm.count("help") || vm.count("task") == 0) {
            if (vm.count("task") == 0 && vm.count("help") == 0)
                std::cout << "MISSING TASK OPTION!" << std::endl;
            std::cerr << cmdlineOptions << std::endl;
            return 0;
        }

        if (vm.count("parameters")) {
            auto par_file = mdtools::open_file(parameters, "expecting parameters file");
            store(parse_config_file(par_file, configFileOptions), vm);
            notify(vm);
        }


        mdtools::create_output_directory(io_options.output_path, io_options.backup);
        mdtools::log_command(io_options.output_path, ac, av);

        std::string log_file = io_options.output_path + "/output.log";

        bool debug = vm.count("debug");
        bool silent = vm.count("silent");
        if (silent)
            mdtools::LOGGER.init(mdtools::log_t::ERROR, log_file);
        else {
            if (debug)
                mdtools::LOGGER.init(mdtools::log_t::DBG, log_file);
            else
                mdtools::LOGGER.init(mdtools::log_t::INFO, log_file);
        }

        mdtools::initialize("Molecular Dynamics Tools.");

        mdtools::show_options(vm);

        io_options.validate();
        switch (mdtools::str2task.at(task)) {

            case mdtools::task_t::PhononDOS :
                phonon_dos.validate();
                mdtools::mainPhononDOS(phonon_dos, io_options,simulation_options);
                break;
            case mdtools::task_t::DynamicStructureFactor :
                dynamic_structure_factor.validate();
                mdtools::mainDynamicStructureFactor(dynamic_structure_factor,io_options,simulation_options);
                break;
            case mdtools::task_t::AxialDistributionHistogram :
                axial_distribution_histogram.validate();
                mdtools::mainAxialDistributionHistogram(axial_distribution_histogram,io_options,simulation_options);
                break;
            default:
                mdtools::LOGGER.error << "Unknown task: " << task << std::endl;
                mdtools::LOGGER.error << "Valid options are: " << std::endl;
                for (auto &item: mdtools::str2task) { mdtools::LOGGER.error << item.first << std::endl; }

        }

        mdtools::LOGGER.close();
        mdtools::finalize(start_time);


    }
    catch (const std::out_of_range &oor) {
        mdtools::LOGGER.error << "Unknown task, valid options are: " << std::endl;
        for (auto &item: mdtools::str2task) { mdtools::LOGGER.error << item.first << std::endl; }
        return return_state;
    }
    catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return return_state;
    }
    catch (...) {
        std::cerr << "Exception of unknown type!" << std::endl;
        return return_state;
    }




    return return_state;

}
