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

#include "definitions.h"
#include "stocc.h"
#include "rhyper.h"
#include <iostream>

int main(int argn, char **argv) {
    StochasticLib1 stocc(0);
    HyperGen gen(0);

    FILE *fp = fopen("deviates", "w+");

    fprintf(fp, "%lld\n", (ULONG) gen.generate(pow(2,16), pow(2,40) - pow(2,16), pow(2,39)));
    fprintf(fp, "%lld\n", (ULONG) gen.generate(pow(2,16), pow(2,60) - pow(2,16), pow(2,59)));
    // for (ULONG i = 16; i < 40; ++i) {
    //     stocc.RandomInit(i);
    //     gen.RandomInit(i);
    //     ULONG N = pow(2,60);
    //     ULONG m = pow(2,59);
    //     ULONG n = pow(2,i);
    //     std::cout << i << std::endl;
    //     fprintf(fp, "%lld\n", (ULONG) gen.generate(n, N - n, m));
    //     // fprintf(fp, "%lld\n", (ULONG) stocc.Hypergeometric(m, n, N));
    // }

    fclose(fp);
}
