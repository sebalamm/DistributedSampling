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
#include <limits>

#include "definitions.h"
#include "randomc.h"
#include "dSFMT.h"

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#ifndef unlikely
#define unlikely(x) __builtin_expect((x),0)
#endif
#ifndef likely
#define likely(x) __builtin_expect((x),1)
#endif

template <typename RandomGenerator = CRandomMersenne, ULONG blocksize = (1 << 24), ULONG dummy = std::numeric_limits<ULONG>::max()>
class HashSampling {
    public:
        HashSampling(ULONG seed, ULONG n) 
            // : gen(seed)
        { 
            // Modification: dSFMT
            dsfmt_init_gen_rand(&dsfmt, seed);
            resizeTable(n);
        }

        void resizeTable(ULONG n) {
            // Table size
            table_lg = 3 + LOG2(n) + isNotPowerOfTwo(n);
            table_size = ipow(2, table_lg);
            hash_table.resize(table_size, dummy);
            
            // Offset for fast indexing
            offset = &(hash_table[0]);
        }

        // See SH subroutine in Ahrens and Dieter
        template <typename F>
        void sample(ULONG N, 
                    ULONG n, 
                    F &&callback) {
            ULONG variate, index, hash_elem;
            ULONG population_lg = (LOG2(N) + isNotPowerOfTwo(N));
            ULONG address_mask = (table_lg >= population_lg) ? 0 : population_lg - table_lg;
            orig_n = n;

            // Modification: dSFMT
            ULONG curr_blocksize = std::max(std::min(n, blocksize), (ULONG)dsfmt_get_min_array_size());
            curr_blocksize += (curr_blocksize & 0x1); // needs to be even
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
                        curr_blocksize += (curr_blocksize & 0x1); // needs to be even
                        // delete[] randblock; randblock = new double[curr_blocksize];
                        dsfmt_fill_array_open_close(&dsfmt, randblock, curr_blocksize);
                        array_index = 0;
                    }
                    variate = N * randblock[array_index++];
                    // Modification: End
                    
                    // variate = N * gen.Random();
                    index = variate >> address_mask; 
                    hash_elem = *(offset + index);    

                    // Table lookup
                    if (likely(hash_elem == dummy)) break; // done
                    else if (hash_elem == variate) continue; // already sampled
                    else {
increment:
                        ++index;
                        index &= (table_size - 1);
                        hash_elem = *(offset + index); 
                        if (hash_elem == dummy) break; // done 
                        else if (hash_elem == variate) continue; // already sampled
                        goto increment; // keep incrementing
                    }
                }
                // Add sample
                *(offset + index) = variate;
                n--;
            }


            // Condense
            ULONG i = 0;
            ULONG j = 0;
            while (i < orig_n) {
                while (*(offset + j) == dummy) j++;
                *(offset + i) = *(offset + j); i++; j++;
            }

            // Exchange sort
            ULONG tmp;   
            for (i = 0; i < orig_n-1; i++) {
                for (j = i+1; j < orig_n; j++) {
                    if (*(offset + i) > *(offset + j)) {
                        tmp = *(offset + i);   
                        *(offset + i) = *(offset + j);
                        *(offset + j) = tmp;
                    }
                }
            }

            // Output in sorted order
            for (i = 0; i < orig_n; i++) callback(*(offset + i) + 1);
            
            clear();
        }

        void clear() {
            // Alternative
            // memset(offset, dummy, sizeof(ULONG)*table_size);
            std::fill(hash_table.begin(), hash_table.begin() + orig_n, dummy);
        }

    private:
        // RandomGenerator gen;
        dsfmt_t dsfmt;

        std::vector<ULONG> hash_table;
        ULONG table_lg, table_size;
        ULONG *offset;
        ULONG orig_n;

        inline ULONG ipow(ULONG base, ULONG exp) {
            ULONG result = 1;
            while (exp) {
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            }
            return result;
        }

        void moveCluster(ULONG index, ULONG variate) {
            ULONG current_elem = *(offset + index);
            ULONG next_index = (index + 1) & (table_size - 1);
            ULONG next_elem = *(offset + next_index);
            while (next_elem != dummy) {
                ++index;
                index &= (table_size - 1);
                // Swap elements
                *(offset + index) = current_elem;
                current_elem = next_elem;
                next_index = (index + 1) & (table_size - 1);
                next_elem = *(offset + next_index);
            } 
            *(offset + next_index) = current_elem;
        }

        inline bool isNotPowerOfTwo(ULONG x) {
            return (x & (x - 1)) != 0;
        }
};

#endif 

