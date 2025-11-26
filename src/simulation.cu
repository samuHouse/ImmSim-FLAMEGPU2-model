#include "simulation.cuh"

SimulationParams parseArguments(int argc, char** argv) {
    SimulationParams params;

    // Maps for parameters of different types
    std::unordered_map<std::string, float*> float_map;
    std::unordered_map<std::string, unsigned int*> uint_map;

    // Float parameters (grid dimensions)
    float_map["-width"]  = &params.width;
    float_map["-height"] = &params.height;

    // Unsigned integer parameters (agent counts, steps)
    uint_map["-B"]     = &params.b_cell_num;
    uint_map["-T"]     = &params.t_cell_num;
    uint_map["-Ag"]    = &params.antigen_num;
    uint_map["-steps"] = &params.steps;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        // Handle float parameters
        if (float_map.count(arg)) {
            if (i + 1 < argc) {
                try {
                    // Convert next argument to float
                    float x = std::stof(argv[i + 1]);
                    if (x > 0) *float_map[arg] = x;
                } catch (...) {
                    std::cerr << "Invalid value for " << arg << "\n";
                }
                i++;
            } else {
                std::cerr << "Missing value after " << arg << "\n";
            }
            continue;
        }

        // Handle unsigned integer parameters
        if (uint_map.count(arg)) {
            if (i + 1 < argc) {
                try {
                    // Convert next argument to unsigned integer
                    unsigned int x = std::stoul(argv[i + 1]);
                    if (x > 0) *uint_map[arg] = x;
                } catch (...) {
                    std::cerr << "Invalid value for " << arg << "\n";
                }
                i++;
            } else {
                std::cerr << "Missing value after " << arg << "\n";
            }
            continue;
        }

        // Handle unknown flag
        std::cerr << "Unknown flag: " << arg << "\n";
    }

    return params;
}