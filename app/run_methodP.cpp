/******************************************************************************
 * run_methodP.cpp 
 *
 * Source of the sampling routine
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <argtable2.h>
#include <mpi.h>

#include <vector>

#include "timer.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "sampling_config.h"
#include "sampling/methodP.h"

int main(int argn, char **argv) {
    // Init MPI
    MPI_Init(&argn, &argv);    
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read command-line args
    SamplingConfig config;
    int ret_code = parse_parameters(argn, argv, 
                                    config); 
    if (ret_code) { MPI_Finalize(); return 0; }

    // Main algorithm
    FILE *fp;
    if (rank == ROOT) {
        std::cout << "sample (n=" << config.n << 
                              ", N=" << config.N << 
                              ", k=" << config.k << 
                              ", s=" << config.seed << 
                              ", p=" << size << ")" << std::endl;
        std::string filename = config.output_file;
        fp = fopen(filename.c_str(), "w+");
    }
    
    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve(config.n);

    // Timers
    timer t;
    t.restart();

    // Compute sample
    ParDivideSampling<> pds(config, config.seed, size);
    pds.sample(config.n,
               1,
               size,
               rank,
               [&](ULONG elem) {
                   // fprintf(fp, "%lld\n", elem);
                   sample.push_back(elem);
               });

    double time_taken = t.elapsed();
    double max_time_taken = 0.0;
    MPI_Reduce(&time_taken, &max_time_taken, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    ULONG samples = sample.size();
    ULONG total_samples = 0;
    MPI_Reduce(&samples, &total_samples, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == ROOT) {
        std::cout << "sampled " << total_samples << " elements" << std::endl;
        std::cout << "total time " << max_time_taken << std::endl;
        fprintf(fp, "total time %f", max_time_taken);
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
