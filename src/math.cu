/****************************************************************************
*
*    HIS Simulator in FLAME GPU
*
*    Copyright (C) 2025  Samuele Casadei
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
#include <time.h>
#include <limits.h>
#include <curand_kernel.h>
#include "math.cuh"
#include "memory.cuh"

__host__ unsigned char randbyte() {
    return rand() % (UCHAR_MAX + 1);
}

__device__ double langevin(double velocity, double force, double mass) {
    return (-LAMBDA * velocity + force) / mass * TIME_FACTOR;
}

__host__ __device__ int hammingdist(unsigned char byte1, unsigned char byte2) {
    unsigned char xored = byte1 ^ byte2;
    // Popcount counts the bits set to 1.
#ifdef __CUDA_ARCH__
    return __popc(xored);
#else
    return __builtin_popcount(xored);
#endif
}
