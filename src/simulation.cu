/****************************************************************************
*
*    ImmSim FLAMEGPU2 model
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

#include "simulation.cuh"

/**
 * This function is used to parse parameters from command line.
 */
SimulationParams parseArguments(int argc, char** argv) {
    SimulationParams params;

    // Two different maps to parse parameters of different types.
    std::unordered_map<std::string, float*> float_map;
    std::unordered_map<std::string, unsigned int*> uint_map;

    // Float parameters (grid dimensions).
    float_map["-width"]  = &params.width;
    float_map["-height"] = &params.height;

    // Unsigned integer parameters (agent counts, steps).
    uint_map["-B"]     = &params.b_cell_num;
    uint_map["-T"]     = &params.t_cell_num;
    uint_map["-Ag"]    = &params.antigen_num;
    uint_map["-steps"] = &params.steps;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];


        if (float_map.count(arg)) {
            if (i + 1 < argc) {
                // Convert next argument to float.
                float x = std::stof(argv[i + 1]);
                if (x > 0) *float_map[arg] = x;
                i++;
            } else {
                std::cerr << "Missing value after " << arg << "\n";
            }
	    // If a float value is found, skip to the next iteration.
            continue;
        }

        if (uint_map.count(arg)) {
            if (i + 1 < argc) {
                // Convert next argument to unsigned integer.
                unsigned int x = std::stoul(argv[i + 1]);
                if (x > 0) *uint_map[arg] = x;
                i++;
            } else {
                std::cerr << "Missing value after " << arg << "\n";
            }
            continue;
        }

        // If an unrecognized flag is set, an error must be thrown.
        std::cerr << "Unknown flag: " << arg << "\n";
    }

    return params;
}