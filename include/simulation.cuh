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
*    Forked and modified by Samuele Casadei, 2025.
*
****************************************************************************/

#include <iostream>
#include <string>
#include <unordered_map>
#include <cstdlib>

#define DEFAULT_WIDTH 500.0f
#define DEFAULT_HEIGHT 500.0f
#define DEFAULT_B_CELL_NUM 200
#define DEFAULT_T_CELL_NUM 200
#define DEFAULT_ANTIGEN_NUM 2000
#define DEFAULT_STEPS 5000

/**
 * Struct containing all the main simulation parameters,
 * initialized at their default value.
 */
struct SimulationParams {
    float width  = DEFAULT_WIDTH;
    float height = DEFAULT_HEIGHT;
    unsigned int b_cell_num = DEFAULT_B_CELL_NUM;
    unsigned int t_cell_num = DEFAULT_T_CELL_NUM;
    unsigned int antigen_num = DEFAULT_ANTIGEN_NUM;
    unsigned int steps = DEFAULT_STEPS;
};

/**
 * This function parses the argv array and overwrites every
 * specified parameter.
 */
SimulationParams parseArguments(int argc, char** argv);
