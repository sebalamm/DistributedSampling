/******************************************************************************
 * run_methodR.cpp 
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

#include <vector>

#include "timer.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "sampling_config.h"
#include "sampling/methodR.h"

int main(int argn, char **argv) {
    // Read command-line args
    SamplingConfig config;
    int ret_code = parse_parameters(argn, argv, 
                                    config); 
    if (ret_code) return 0;

    // Main algorithm
    FILE *fp;
    std::cout << "sample (n=" << config.n << 
                          ", N=" << config.N << 
                          ", k=" << config.k << 
                          ", s=" << config.seed << ")" << std::endl;
    std::string filename = config.output_file;
    fp = fopen(filename.c_str(), "w+");

    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve(config.n);

    // Timers
    timer t;
    t.restart();

    // Compute sample
    HashSampling<> hs(config.seed);
    hs.resizeTable(config.N / config.n * config.k, config.k);
    SeqDivideSampling<> sds(hs, config.k, config.seed);
    sds.sample(config.N,
               config.n,
               [&](ULONG elem) {
                   // fprintf(fp, "%lld\n", elem);
                   sample.push_back(elem);
               });

    std::cout << "sampled " << sample.size() << " elements" << std::endl;
    std::cout << "total time " << t.elapsed() << std::endl;
    fprintf(fp, "total time %f", t.elapsed());
    fclose(fp);
}

