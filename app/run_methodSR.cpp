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

#include <vector>
#include <algorithm>

#include "timer.h"
#include "macros_assertions.h"
#include "sampling_config.h"
#include "sampling/methodR.h"
#include "tools/benchmark.h"
#include "tools/arg_parser.h"

int main(int argn, char **argv) {
    // Read command-line args
    SamplingConfig config;

    arg_parser args(argn, argv);
    config.N = 1 << args.get<ULONG>("N", 30);
    config.n = 1 << args.get<ULONG>("n", 20); // sample size
    config.k = 1 << args.get<ULONG>("k", 10); // base case
    config.seed = args.get<ULONG>("seed", 1);
    config.iterations = args.get<ULONG>("i", 1);
    config.output_file = args.get<std::string>("output", "tmp");

    // Main algorithm
    FILE *fp;
    std::cout << "sample (n=" << config.n << 
                          ", N=" << config.N << 
                          ", k=" << config.k << 
                          ", s=" << config.seed <<
                          ", i=" << config.iterations << ")" << std::endl;
    std::string filename = config.output_file;
    fp = fopen(filename.c_str(), "w+");

    // Resulting samples
    std::vector<ULONG> sample;
    sample.reserve(config.n);

    // Statistics
    timer t;
    statistics stats;

    std::cout << "warmup" << std::endl;
    for (ULONG iteration = 0; iteration < std::min((ULONG)100, config.iterations); ++iteration) {
        sample.clear();

        // Compute sample
        SortedHashSampling<> hs(config.seed + iteration, config.k);
        SeqDivideSampling<StochasticLib1, SortedHashSampling<>> sds(hs, config.k, config.seed + iteration);
        sds.sample(config.N,
                   config.n,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                   });

        
        if (!std::is_sorted(sample.begin(), sample.end())) std::cout << "not sorted!" << std::endl;
        if (sample.size() != config.n) std::cout << "wrong size " << sample.size() << "!" << std::endl;
    }

    std::cout << "measurements" << std::endl;
    for (ULONG iteration = 0; iteration < config.iterations; ++iteration) {
        sample.clear();
        t.restart();

        // Compute sample
        SortedHashSampling<> hs(config.seed + iteration, config.k);
        SeqDivideSampling<StochasticLib1, SortedHashSampling<>> sds(hs, config.k, config.seed + iteration);
        sds.sample(config.N,
                   config.n,
                   [&](ULONG elem) {
                       // fprintf(fp, "%lld\n", elem);
                       sample.push_back(elem);
                   });

        double time = t.elapsed();
        stats.push(time);    
    }

    std::cout << "RESULT runner=SR" 
              << " time=" << stats.avg()
              << " stddev=" << stats.stddev()
              << " iterations=" << config.iterations << std::endl;
    fprintf(fp, "RESULT runner=SR time=%f stddev=%f iterations=%llu\n", stats.avg(), stats.stddev(), config.iterations);
    fclose(fp);
}

