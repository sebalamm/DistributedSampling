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
#include <iostream>

int main(int argn, char **argv) {
    StochasticLib1 stocc(0);

    for (ULONG i = 0; i < 100; ++i) {
        stocc.RandomInit(i);
        ULONG N = pow(2,60);
        // ULONG N = pow(2,60);
        ULONG m = N/2;
        ULONG n = pow(2,40);
        // ULONG n = pow(2,10);
        std::cout << stocc.Hypergeometric(n, m, N) << std::endl;
    }
}
