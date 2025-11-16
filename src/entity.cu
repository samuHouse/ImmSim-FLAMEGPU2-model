/****************************************************************************
*    
*    HIS Simulator in C/Cuda C++
*
*    Copyright (C) 2025  Daniel Pellanda
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <https://www.gnu.org/licenses/>
*
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cmath>
#include <climits>
#include <curand_kernel.h>
#include "math.cuh"
#include "entity.cuh"
#include "memory.cuh"

__device__ double affinity_potential(unsigned char receptor1, unsigned char receptor2) {
    int dist = hammingdist(receptor1, receptor2);
    if (dist < AFFINITY_MIN)
        return 0.0;
    return std::pow(BIND_CHANCE, (double)(dist - CHAR_BIT) / (double)(AFFINITY_MIN - CHAR_BIT));
}

__device__ bool can_entities_bind(const unsigned char* receptor1, const unsigned char* receptor2) {
    for (int i = 0; i < RECEPTOR_SIZE; i++) {
        /* Only one pair of receptors needs to bind. */
        if (randdouble() < affinity_potential(receptor1[i], receptor2[i]))
            return true;
    }
    return false;
}

__host__ __device__ void hypermutation(unsigned char* receptor) {
    #ifdef POISSON_MUTATION
        // Poisson distribution
    #ifdef __CUDA_ARCH__
        curandState state;
        curand_init(clock64(), threadIdx.x + blockIdx.x * blockDim.x, 0, &state);
    #endif
        /* Choose how many and which bits are going to get changed
           and store the indexes in an array. */
        int num_bits = rand() % (CHAR_BIT * RECEPTOR_SIZE); // extract bits to change
        int *positions = (int*)malloc(num_bits * sizeof(int));
        for (int i = 0; i < num_bits; i++) {
    #ifdef __CUDA_ARCH__
	    positions[i] = curand(&state) % (CHAR_BIT * RECEPTOR_SIZE);
    #else
            positions[i] = rand() % (CHAR_BIT * RECEPTOR_SIZE); // extract index to change
    #endif
        }

        /* Invert the value of every bit in the position of each index contained in the array. */
        for (int i = 0; i < num_bits; i++) {
            int pos = positions[i] % CHAR_BIT;
            int index = positions[i] / CHAR_BIT;
            // check the value of a given bit in a byte knowing its position.
            bool set = (receptor[index] & (1 << pos)) != 0;
            if (set) {
                // when the bit is 1, force it to zero through a bitwise AND.
                receptor[index] &= ~(1 << pos);
            } else {
                // when the bit is 0, force it to one through a bitwise OR.
                receptor[index] |= (1 << pos);
            }
        }
        free(positions);
    #else
        // Binomial distribution
        for (int i = 0; i < RECEPTOR_SIZE; i++) {
            for (int j = 0; j < BITS_IN_A_BYTE; j++) {
                double random = randdouble();
                if (random < MUTATION_CHANCE) {
                    bool set = !getbit(entity->receptor[i], j);
                    setbit(&entity->receptor[i], set, j);
                }
            }
        }
    #endif
}
