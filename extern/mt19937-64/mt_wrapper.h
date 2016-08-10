/******************************************************************************
 * mt_wrapper.h
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


#ifndef _MT_WRAPPER_H_
#define _MT_WRAPPER_H_

extern "C" {
    #include "mt64.h"
    void init_genrand64(unsigned long long seed);
    void init_by_array64(unsigned long long init_key[], 
                 unsigned long long key_length);
    unsigned long long genrand64_int64(void);
    double genrand64_real2(void);
}

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(const char *ErrorText);// System-specific error reporting (userintf.cpp)

class MTWrapper {
    public:
        MTWrapper() {};

        MTWrapper(unsigned long long seed) {
            RandomInit(seed);
        };

        void RandomInit(unsigned long long seed) {
            init_genrand64(seed);
        }

        void RandomInitByArray(unsigned long long seeds[], unsigned long long NumSeeds) {
            init_by_array64(seeds, NumSeeds);
        }

        unsigned long long BRandom() {
            return genrand64_int64();
        }

        double Random() {
            return genrand64_real2();
        }

        unsigned long long IRandom(unsigned long long min, unsigned long long max) {
            if (max == min) return min;
            unsigned long long r = (unsigned long long)((double)(unsigned long long)(max - min + 1) * Random() + min); 
            if (r > max) r = max;
            return r;
        }
};

#endif
