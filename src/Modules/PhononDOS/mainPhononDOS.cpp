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


#include <fftw3.h>
#include "mainPhononDOS.h"
#include "trajectoryReader.h"
#include "logger.h"

namespace mdtools {


    void mainPhononDOS(phonon_dos_options_t phonon_dos, const io_options_t &io_options,
                       const simulation_options_t &simulation_options) {

        LOGGER.info << "main Phonon DOS" << std::endl;

        auto trajectory = trajectoryReader(io_options.trajectory_input_file,io_options.coordinates_input_file).get(simulation_options.time_step,
                                                                                 simulation_options.start_iteration,
                                                                                 simulation_options.delta_iteration,
                                                                                 simulation_options.end_iteration);

        std::valarray<double> vaf = std::valarray<double>(0.0, trajectory[0].velocity_x.size());
        std::valarray<double> norm = std::valarray<double>(0.0, trajectory[0].velocity_x.size());

        auto n = trajectory[0].velocity_x.size();

        auto inorm = 1.0 / trajectory.size();

        for (auto &item: trajectory) {
            for (int i = 0; i < n; i++) {
                auto vxi = item.velocity_x[i];
                auto vyi = item.velocity_y[i];
                auto vzi = item.velocity_z[i];
                auto in = 1.0 / (n - i);
                for (int j = i; j < n; j++) {
                    vaf[j] += (vxi * item.velocity_x[j] + vyi * item.velocity_y[j] + vzi * item.velocity_z[j]) * inorm *
                              in;
                    norm[j] = n - i;
                }
            }
        }


        std::ofstream file;
        file.open("vacf.csv", std::ios::out);
        LOGGER.debug << "writing output" << std::endl;
        file << "time,my_vacf,norm" << std::endl;
        for (int i = 0; i < vaf.size(); i++) {
            file << i * 0.1 << "," << vaf[i] << "," << norm[i] << std::endl;
        }

        file.close();

        /*
        fftw_complex *in, *out;
        fftw_plan p;

        auto N=100;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        for(int i=0; i < N; i++){
            in[0][i]=sin(3.141593*i/4);
        }
        fftw_execute(p);

        for(int i=0; i < N; i++){
            LOGGER.debug << sqrt(out[0][i]*out[0][i] + out[1][i]*out[1][i]) << std::endl;
        }

        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);
       */

    }

}