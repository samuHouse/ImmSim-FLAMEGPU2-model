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

struct SimulationParams {
    float width  = DEFAULT_WIDTH;
    float height = DEFAULT_HEIGHT;
    unsigned int b_cell_num = DEFAULT_B_CELL_NUM;
    unsigned int t_cell_num = DEFAULT_T_CELL_NUM;
    unsigned int antigen_num = DEFAULT_ANTIGEN_NUM;
    unsigned int steps = DEFAULT_STEPS;
};

SimulationParams parseArguments(int argc, char** argv);
