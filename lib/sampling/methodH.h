/******************************************************************************
 * methodH.h
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


#ifndef _METHOD_H_H_
#define _METHOD_H_H_

#include <vector>
#include "definitions.h"
#include "randomc.h"
#include "dSFMT.h"

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define unlikely(x) __builtin_expect((x),0)
#define likely(x) __builtin_expect((x),1)

template <typename RandomGenerator = CRandomMersenne, ULONG blocksize = (1 << 24)>
class HashSampling {
    public:
        HashSampling(ULONG seed) 
            // : gen(seed)
        { 
            // Modification: dSFMT
            dsfmt_init_gen_rand(&dsfmt, seed);
        }

        void resizeTable(ULONG N, ULONG n) {
            ULONG table_lg = 3 + LOG2(n);
            address_mask = LOG2(N) - table_lg;

            // Table size
            table_size = ipow(2, table_lg);
            hash_table.resize(table_size, 0);
            indices.reserve(table_size);
            
            // Offset for fast indexing
            offset = &(hash_table[0]);
        }

        // See SH subroutine in Ahrens and Dieter
        template <typename F>
        void sample(ULONG N, 
                    ULONG n, 
                    F &&callback) {
            ULONG variate, index, hash_elem;

            // Modification: dSFMT
            ULONG curr_blocksize = std::max(std::min(n, blocksize), (ULONG)dsfmt_get_min_array_size());
            double *randblock = new double[curr_blocksize];
            dsfmt_fill_array_open_close(&dsfmt, randblock, curr_blocksize);
            ULONG array_index = 0;
            // Modification: End

            while (n > 0) {
                while (true) {
                    // Take sample
  
                    // Modification: dSFMT
                    if (array_index >= curr_blocksize) {
                        curr_blocksize = std::max(std::min(n, blocksize), (ULONG)dsfmt_get_min_array_size());
                        delete[] randblock; randblock = new double[curr_blocksize];
                        dsfmt_fill_array_open_close(&dsfmt, randblock, curr_blocksize);
                        array_index = 0;
                    }
                    variate = N * randblock[array_index++] + 1;
                    // Modification: End
                    
                    // variate = N * gen.Random() + 1;
                    index = variate >> address_mask; 
                    hash_elem = *(offset + index);    

                    // Table lookup
                    if (likely(hash_elem == 0)) break; // done
                    else if (hash_elem == variate) continue; // already sampled
                    else {
increment:
                        ++index;
                        if (unlikely(index >= table_size)) index = 0; // restart probing
                        hash_elem = *(offset + index); 
                        if (hash_elem == 0) break; // done 
                        else if (hash_elem == variate) continue; // already sampled
                        goto increment; // keep incrementing
                    }
                }
                // Add sample
                *(offset + index) = variate;
                indices.push_back(index);
                callback(variate);
                n--;
            }

            clear();
        }

        bool isEmpty() {
            return indices.empty();
        }

        void clear() {
            for (ULONG index : indices) hash_table[index] = 0; 
            indices.clear();
            // Alternative
            // memset(offset, 0, sizeof(ULONG)*table_size);
            // std::fill(hash_table.begin(), hash_table.end(), 0);
        }

    private:
        // RandomGenerator gen;
        dsfmt_t dsfmt;

        std::vector<ULONG> indices;
        std::vector<ULONG> hash_table;
        ULONG table_size;
        ULONG address_mask;
        ULONG *offset;

        ULONG ipow(ULONG base, ULONG exp) {
            ULONG result = 1;
            while (exp) {
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            }
            return result;
        }
};

#endif 

