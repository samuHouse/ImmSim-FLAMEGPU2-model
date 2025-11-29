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

#include <string>
#include <random>
#include <chrono>
#include "math.cuh"
#include "entity.cuh"
#include "simulation.cuh"
#include "flamegpu/flamegpu.h"

#ifdef __PLOT__

// Not needed python libraries in matplotlib.
#define WITHOUT_NUMPY
#define WITHOUT_OPENCV
#include "matplotlibcpp.h"
// necessary in order to use plot functions.
namespace plt = matplotlibcpp;

#endif

/************************************ AGENT FUNCTION DEFINITION *************************************
 *
 *	Agent functions run on the GPU, every agent runs it in an isolated context, each agent
 *	can only access its local data. In order to exchange information with other agents,
 *	each agent can make use of input/output messages defined in the signature of each function.
 *	Agent functions are executed for each agent ( they're added to) at each timestep of the
 *	simulation.
 */

FLAMEGPU_AGENT_FUNCTION(antigen_send_pos, flamegpu::MessageNone, flamegpu::MessageSpatial2D) {
    // Get the agent's position
    const float x = FLAMEGPU->getVariable<float>("x");
    const float y = FLAMEGPU->getVariable<float>("y");
    unsigned char receptor[RECEPTOR_SIZE];
    for (int i=0; i<RECEPTOR_SIZE; i++) {
        receptor[i] = FLAMEGPU->getVariable<unsigned char, RECEPTOR_SIZE>("epitope", i);
        FLAMEGPU->message_out.setVariable<unsigned char, RECEPTOR_SIZE>("receptor", i, receptor[i]);
    }
    // Send position message
    FLAMEGPU->message_out.setVariable<float>("x", x);
    FLAMEGPU->message_out.setVariable<float>("y", y);
    FLAMEGPU->message_out.setVariable<unsigned int>("id", FLAMEGPU->getID());

    // Behaviour goes here
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(antibody_interaction, flamegpu::MessageSpatial2D, flamegpu::MessageSpatial2D) {
    uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    if (!has_interacted) {
    const float RADIUS = FLAMEGPU->message_in.radius();
	const float x = FLAMEGPU->getVariable<float>("x");
	const float y = FLAMEGPU->getVariable<float>("y");
        for (auto entity : FLAMEGPU->message_in(x, y)) {
            // if message is in range and not out of bounds
            const float x2 = entity.getVariable<float>("x");
            const float y2 = entity.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // check if antigen is in range for interaction.
                float x21 = x2 - x;
                float y21 = y2 - y;
                const float separation = sqrt( x21*x21 + y21*y21 );
                if (separation < RADIUS && separation > 0.0f) {
                    // Antibody binds to antigen, sending kill message.
                    FLAMEGPU->message_out.setVariable<float>("x", entity.getVariable<float>("x"));
                    FLAMEGPU->message_out.setVariable<float>("y", entity.getVariable<float>("y"));
                    FLAMEGPU->message_out.setVariable<unsigned int>("id", entity.getVariable<unsigned int>("id"));
                    // Mark as interacted to prevent multiple bindings.
                    FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
                    break; // Exit after first successful interaction.
                }
            }
        }
    }
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(antigen_receive_kill, flamegpu::MessageSpatial2D, flamegpu::MessageNone) {
    const float RADIUS = FLAMEGPU->message_in.radius();
	const float x = FLAMEGPU->getVariable<float>("x");
	const float y = FLAMEGPU->getVariable<float>("y");
    for (auto msg : FLAMEGPU->message_in(x, y)) {
        // If message is in range and not out of bounds
        const float x2 = msg.getVariable<float>("x");
        const float y2 = msg.getVariable<float>("y");
        if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
            0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
        )
        {
            // Check if antigen is in range for interaction.
            float x21 = x2 - x;
            float y21 = y2 - y;
            const float separation = sqrt( x21*x21 + y21*y21 );
            if (separation < RADIUS && separation > 0.0f) {
                if (msg.getVariable<flamegpu::id_t>("id") == FLAMEGPU->getID()) {
                    return flamegpu::DEAD; // Antigen is killed
                }
            }
        }
    }
    return flamegpu::ALIVE; // Antigen remains alive
}

FLAMEGPU_AGENT_FUNCTION(b_cell_binding, flamegpu::MessageSpatial2D, flamegpu::MessageSpatial2D) {
    uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    CellState state = static_cast<CellState>(FLAMEGPU->getVariable<int>("cell_state"));
    if (!has_interacted && state == CS_INTERNALIZED) {
        const float RADIUS = FLAMEGPU->message_in.radius();
	const float x = FLAMEGPU->getVariable<float>("x");
	const float y = FLAMEGPU->getVariable<float>("y");
	const flamegpu::AgentRandom &rng = FLAMEGPU->random;
        for (auto entity : FLAMEGPU->message_in(x, y)) {
            // If message is in range and not out of bounds
            const float x2 = entity.getVariable<float>("x");
            const float y2 = entity.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // Check if antigen is in range for interaction.
                float x21 = x2 - x;
                float y21 = y2 - y;
                const float separation = sqrt( x21*x21 + y21*y21 );
                if (separation < RADIUS && separation > 0.0f) {
                    unsigned char ag_receptor[RECEPTOR_SIZE];
                    unsigned char receptor[RECEPTOR_SIZE];
		    for (int i=0; i<RECEPTOR_SIZE; i++) {
			ag_receptor[i] = entity.getVariable<unsigned char, RECEPTOR_SIZE>("receptor", i);
			receptor[i] = FLAMEGPU->getVariable<unsigned char, RECEPTOR_SIZE>("receptor", i);

		    }
                    if (can_entities_bind(receptor, ag_receptor, rng)) {
                        FLAMEGPU->setVariable<int>("cell_state", CS_ACTIVE);
                        FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
                        FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
                        FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
                        FLAMEGPU->message_out.setVariable<unsigned int>("id", entity.getVariable<unsigned int>("id"));
                        break; // Exit after first successful binding
                    }
                }
            }
        }
    }
    // Interaction logic goes here
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(b_cell_send_pos, flamegpu::MessageNone, flamegpu::MessageSpatial2D) {
    // Get the agent's position
    const float x = FLAMEGPU->getVariable<float>("x");
    const float y = FLAMEGPU->getVariable<float>("y");
    const unsigned int id = FLAMEGPU->getID();
    const int state = static_cast<CellState>(FLAMEGPU->getVariable<int>("cell_state"));
    const uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    // Send position message
    FLAMEGPU->message_out.setVariable<float>("x", x);
    FLAMEGPU->message_out.setVariable<float>("y", y);
    FLAMEGPU->message_out.setVariable<unsigned int>("id", id);
    FLAMEGPU->message_out.setVariable<int>("cell_state", state);
    FLAMEGPU->message_out.setVariable<uint8_t>("has_interacted", has_interacted);
    // Behaviour goes here
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(t_cell_send_stimulus, flamegpu::MessageSpatial2D, flamegpu::MessageSpatial2D) {
    const uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    if (!has_interacted) {
        const float RADIUS = FLAMEGPU->message_in.radius();
        const float x = FLAMEGPU->getVariable<float>("x");
        const float y = FLAMEGPU->getVariable<float>("y");
        for (auto& b_cell : FLAMEGPU->message_in(x, y)) {
            // If message is in range and not out of bounds
            const float x2 = b_cell.getVariable<float>("x");
            const float y2 = b_cell.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // Check if antigen is in range for interaction.
                float x21 = x2 - x;
                float y21 = y2 - y;
                const float separation = sqrt( x21*x21 + y21*y21 );
                if (separation < RADIUS && separation > 0.0f) {
                    const int state = b_cell.getVariable<int>("cell_state");
                    const uint8_t has_interacted = b_cell.getVariable<uint8_t>("has_interacted");
                    // If any nearby B-cell is in the active state and hasn't interacted yet...
                    if (!has_interacted && state == CS_ACTIVE) {
                        const float x = FLAMEGPU->getVariable<float>("x");
                        const float y = FLAMEGPU->getVariable<float>("y");
                        const unsigned int id = b_cell.getVariable<unsigned int>("id");
                        FLAMEGPU->message_out.setVariable<float>("x", x);
                        FLAMEGPU->message_out.setVariable<float>("y", y);
                        FLAMEGPU->message_out.setVariable<unsigned int>("id", id);
                        // Update interaction flag
                        FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
                        break;
                    }
                }
            }
        }
    }
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(b_cell_receive_stimulus, flamegpu::MessageSpatial2D, flamegpu::MessageNone) {
    const uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    const CellState state = static_cast<CellState>(FLAMEGPU->getVariable<int>("cell_state"));
    if (!has_interacted && state == CS_ACTIVE) {
        const float RADIUS = FLAMEGPU->message_in.radius();
        const float x = FLAMEGPU->getVariable<float>("x");
        const float y = FLAMEGPU->getVariable<float>("y");
        for (auto stimulus : FLAMEGPU->message_in(x, y)) {
            // If message is in range and not out of bounds
            const float x2 = stimulus.getVariable<float>("x");
            const float y2 = stimulus.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // Check if antigen is in range for interaction.
                float x21 = x2 - x;
                float y21 = y2 - y;
                const float separation = sqrt( x21*x21 + y21*y21 );
                if (separation < RADIUS && separation > 0.0f) {
                    if (stimulus.getVariable<flamegpu::id_t>("id") == FLAMEGPU->getID()) {
                        FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
                        FLAMEGPU->setVariable<int>("cell_state", CS_STIMULATED);
                        break; // Exit after first stimulus received
                    }
                }
            }
        }
    }
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(spawn_antibodies, flamegpu::MessageNone, flamegpu::MessageNone) {
    const CellState state = static_cast<CellState>(FLAMEGPU->getVariable<int>("cell_state"));
    const uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    const float range = 5,
	min_range = 0.5; // Define a range for spawning antibodies around the B-cell
    int spawned_count = 0;
    if (state == CS_STIMULATED && !has_interacted) {
        const float x = FLAMEGPU->getVariable<float>("x");
        const float y = FLAMEGPU->getVariable<float>("y");
	const flamegpu::AgentRandom &rng = FLAMEGPU->random;
	
	for (int i=0; i < AB_CREATED_PER_CELL; i++) {
	    // Random offset from the parent B-cell.
	    float dx = FLAMEGPU->random.uniform<float>(min_range, range);
	    float dy = FLAMEGPU->random.uniform<float>(min_range, range);
	    // Randomly decide if the offset is negative.
	    dx *= (FLAMEGPU->random.uniform<float>() < 0.5f) ? -1 : 1;
	    dy *= (FLAMEGPU->random.uniform<float>() < 0.5f) ? -1 : 1;
	    // Antibody spawns.
	    FLAMEGPU->agent_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x") + dx);
	    FLAMEGPU->agent_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y") + dy);

	}
        FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
        FLAMEGPU->setVariable<int>("cell_state", CS_INTERNALIZED);

	// Mutate the receptor.
	unsigned char receptor[RECEPTOR_SIZE];
	for (int i=0; i<RECEPTOR_SIZE; i++) receptor[i] = FLAMEGPU->getVariable<unsigned char, RECEPTOR_SIZE>("receptor", i);
	hypermutation(receptor, rng);
	for (int i=0; i<RECEPTOR_SIZE; i++) FLAMEGPU->setVariable<unsigned char, RECEPTOR_SIZE>("receptor", i, receptor[i]);
    }
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(diffuse_entities, flamegpu::MessageNone, flamegpu::MessageNone) {
    // Load grid dimensions off environment.
    const float width  = FLAMEGPU->environment.getProperty<float>("width");
    const float height = FLAMEGPU->environment.getProperty<float>("height");

    // Fetch the current agent information needed for langevin.
    float x = FLAMEGPU->getVariable<float>("x");
    float y = FLAMEGPU->getVariable<float>("y");
    float vx = FLAMEGPU->getVariable<float>("vx");
    float vy = FLAMEGPU->getVariable<float>("vy");
    const float mass = FLAMEGPU->getVariable<float>("mass");
    const float EPS = 1e-3f;

    // Box-Muller
    double rx = FLAMEGPU->random.normal<double>(); // mean 0, std 1
    double ry = FLAMEGPU->random.normal<double>();

    // Langevin
    vx += langevin(vx, rx, mass);
    vy += langevin(vy, ry, mass);

    // Proposed new position.
    float new_x = x + vx * TIME_FACTOR;
    float new_y = y + vy * TIME_FACTOR;

    // Reflection on bounds.
    if (new_x < 0.0f) {
    	new_x = -new_x;
    	vx = -vx;
    } else if (new_x > width) {
    	new_x = 2*width - new_x;
    	vx = -vx;
    }

    if (new_y < 0.0f) {
    	new_y = -new_y;
    	vy = -vy;
    } else if (new_y > height) {
    	new_y = 2*height - new_y;
   	vy = -vy;
    }

    // Update positions.
    FLAMEGPU->setVariable<float>("x", new_x);
    FLAMEGPU->setVariable<float>("y", new_y);

    // Update velocities.
    FLAMEGPU->setVariable<float>("vx", vx);
    FLAMEGPU->setVariable<float>("vy", vy);

    // Reset interaction flag.
    EntityType type = static_cast<EntityType>(FLAMEGPU->getVariable<int>("type"));
    if (type != AG_MOLECOLE) {
        FLAMEGPU->setVariable<uint8_t>("has_interacted", false);
    }

    return flamegpu::ALIVE;
}

/************************************ INIT FUNCTION DEFINITION *************************************
 *
 *	Each simulation has an init function, used to populate the environment with agents,
 *	this type of function runs only once, at the beginning of the simulation.
 *	Runs on CPU (host).
 */

FLAMEGPU_INIT_FUNCTION(init_agents) {
    const unsigned int width = FLAMEGPU->environment.getProperty<float>("width");
    const unsigned int height = FLAMEGPU->environment.getProperty<float>("height");
    const unsigned int b_cell_num = FLAMEGPU->environment.getProperty<unsigned int>("B_cell_num");
    const unsigned int t_cell_num = FLAMEGPU->environment.getProperty<unsigned int>("T_cell_num");
    const unsigned int antigen_num = FLAMEGPU->environment.getProperty<unsigned int>("antigen_num");
    const unsigned int total_agents = b_cell_num + t_cell_num + antigen_num;

    // Create rng to shuffle the types array.
    std::random_device rd;
    std::mt19937 gen(rd());

    // Shuffle types array to populate the grid.
    std::vector<EntityType> types;
    types.reserve(total_agents);

    types.insert(types.end(), b_cell_num, B_CELL);
    types.insert(types.end(), t_cell_num, T_CELL);
    types.insert(types.end(), antigen_num, AG_MOLECOLE);

    std::shuffle(types.begin(), types.end(), gen);

    auto b_cell = FLAMEGPU->agent("B-Lymphocyte");
    auto t_cell = FLAMEGPU->agent("T-Lymphocyte");
    auto antigen = FLAMEGPU->agent("Antigen");
    auto antibody = FLAMEGPU->agent("Antibody");

    unsigned int counter = 0;

    // How many agents in the grid
    unsigned int cols = static_cast<unsigned int>(std::ceil(std::sqrt(total_agents)));
    unsigned int rows = static_cast<unsigned int>(std::ceil(static_cast<float>(total_agents) / cols));

    float cell_width = static_cast<float>(width) / cols;
    float cell_height = static_cast<float>(height) / rows;

    // RNG per posizionamento
    std::uniform_real_distribution<float> offset(0.0f, 1.0f);

    for (unsigned int i = 0; i < rows && counter < total_agents; ++i) {
        for (unsigned int j = 0; j < cols && counter < total_agents; ++j) {
	    // Spread agents equally with a random offset.
	    float x = j * cell_width  + offset(gen) * cell_width;
            float y = i * cell_height + offset(gen) * cell_height;

	    EntityType type = types[counter++];
            switch(type) {
                case B_CELL: {
                    auto instance = b_cell.newAgent();
                    //b++;
                    instance.setVariable<int>("cell_state", CS_INTERNALIZED);
                    instance.setVariable<int>("type", B_CELL);
                    instance.setVariable<float>("mass", 0.2f);
                    instance.setVariable<uint8_t>("has_interacted", false);
                    for (unsigned int k = 0; k < RECEPTOR_SIZE; ++k)
                        instance.setVariable<unsigned char, RECEPTOR_SIZE>("receptor", k, randbyte());
                    instance.setVariable<float>("x", x);
                    instance.setVariable<float>("y", y);
                    instance.setVariable<float>("vx", 0.0f);
                    instance.setVariable<float>("vy", 0.0f);
                    break;
                }

                case T_CELL: {
                    auto instance = t_cell.newAgent();
                    //t++;
                    instance.setVariable<int>("cell_state", CS_ACTIVE);
                    instance.setVariable<int>("type", T_CELL);
                    instance.setVariable<float>("mass", 0.2f);
                    instance.setVariable<uint8_t>("has_interacted", false);
                    instance.setVariable<float>("x", x);
                    instance.setVariable<float>("y", y);
                    instance.setVariable<float>("vx", 0.0f);
                    instance.setVariable<float>("vy", 0.0f);
                    break;
                }

                case AG_MOLECOLE: {
                    auto instance = antigen.newAgent();
                    //ag++;
                    instance.setVariable<int>("cell_state", CS_ACTIVE);
                    instance.setVariable<int>("type", AG_MOLECOLE);
                    instance.setVariable<float>("mass", 0.1f);
                    for (unsigned int k = 0; k < RECEPTOR_SIZE; ++k)
                        instance.setVariable<unsigned char, RECEPTOR_SIZE>("epitope", k, randbyte());
                    instance.setVariable<float>("x", x);
                    instance.setVariable<float>("y", y);
                    instance.setVariable<float>("vx", 0.0f);
                    instance.setVariable<float>("vy", 0.0f);
                    break;
                }

                default: break;
            }
        }
    }
}

/************************************ HOST FUNCTION DEFINITION *************************************
 *
 *	Host functions are designed to retrieve/edit the simulation data.
 *	This type of function is slower compared to agent functions, and,
 *	usually, isn't supposed to run at every timestep.
 *	Runs on CPU (host).
 */

#ifdef __PLOT__

FLAMEGPU_HOST_FUNCTION(plot_agents) {
    const unsigned int step = FLAMEGPU->getStepCounter();
    const unsigned int interval = 500;
    const char* home = std::getenv("HOME");
    std::string plots_path = std::string(home) + "/TesiImmSim/plots/";

    if (step % interval == 0) {
        std::vector<float> x_b, y_b;
        std::vector<float> x_t, y_t;
        std::vector<float> x_ag, y_ag;
        std::vector<float> x_ab, y_ab;

        // Get agents' positions through population data.
        flamegpu::DeviceAgentVector b_cells = FLAMEGPU->agent("B-Lymphocyte").getPopulationData();
	flamegpu::DeviceAgentVector t_cells = FLAMEGPU->agent("T-Lymphocyte").getPopulationData();
	flamegpu::DeviceAgentVector antigens = FLAMEGPU->agent("Antigen").getPopulationData();
	flamegpu::DeviceAgentVector antibodies  = FLAMEGPU->agent("Antibody").getPopulationData();

        // Collect positions.
        for (auto b : b_cells) {
            // Append positions to respective vectors.
            x_b.push_back(b.getVariable<float>("x"));
            y_b.push_back(b.getVariable<float>("y"));
        }
        for (auto t : t_cells) {
            x_t.push_back(t.getVariable<float>("x"));
            y_t.push_back(t.getVariable<float>("y"));
        }
        for (auto ag : antigens) {
            x_ag.push_back(ag.getVariable<float>("x"));
            y_ag.push_back(ag.getVariable<float>("y"));
        }
        for (auto ab : antibodies) {
            x_ab.push_back(ab.getVariable<float>("x"));
            y_ab.push_back(ab.getVariable<float>("y"));
        }

        // Plotting.
        plt::clf();
        plt::scatter(x_b, y_b, 10.0, {{"color", "orchid"}, {"label", "B-Cells"}});
        plt::scatter(x_t, y_t, 10.0, {{"color", "green"}, {"label", "T-Cells"}});
        plt::scatter(x_ag, y_ag, 10.0, {{"color", "red"}, {"label", "Antigens"}});
        plt::scatter(x_ab, y_ab, 10.0, {{"color", "blue"}, {"label", "Antibodies"}});
        plt::title("Step " + std::to_string(step));
        plt::xlim(0, static_cast<int>(FLAMEGPU->environment.getProperty<float>("width")));
        plt::ylim(0, static_cast<int>(FLAMEGPU->environment.getProperty<float>("height")));
        plt::save( plots_path + "/step_" + std::to_string(step) + ".png" );
    }
}

#endif

FLAMEGPU_HOST_FUNCTION(reinsert_antigens) {

    unsigned int total_steps = FLAMEGPU->environment.getProperty<unsigned int>("total_steps"),
    step = FLAMEGPU->getStepCounter();
    if (step == (total_steps / 2)) {
    	// Gather environment variables.
    	const unsigned int width  = (unsigned int)FLAMEGPU->environment.getProperty<float>("width");
    	const unsigned int height = (unsigned int)FLAMEGPU->environment.getProperty<float>("height");
    	const unsigned int antigen_num = FLAMEGPU->environment.getProperty<unsigned int>("antigen_num");

    	// Fetch agent population in order to mark occupied cells.
    	unsigned int bcells = FLAMEGPU->agent("B-Lymphocyte").count();
    	unsigned int tcells = FLAMEGPU->agent("T-Lymphocyte").count();
    	unsigned int abs = FLAMEGPU->agent("Antibody").count();
    	unsigned int ags = FLAMEGPU->agent("Antigen").count();

    	// Also counting reinserted antigens.
    	const unsigned int total_agents = bcells + tcells + ags + abs + antigen_num;

    	unsigned int cols = static_cast<unsigned int>(std::ceil(std::sqrt(total_agents)));
    	unsigned int rows = static_cast<unsigned int>(std::ceil((float)total_agents / cols));
    	float cell_width = (float)width / cols,
    	cell_height = (float)height / rows;

    	// Gather antigen proxy, to create more.
    	auto antigen = FLAMEGPU->agent("Antigen");

    	// Shuffle free_positions.
    	std::random_device rd;
    	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> distX(0.0f, width);
	std::uniform_real_distribution<float> distY(0.0f, height);

	for (unsigned int k = 0; k < antigen_num; ++k) {
    	    float x, y;
            // posizione completamente casuale
            x = distX(gen);
            y = distY(gen);

    	    // spawn dell’agente
    	    auto instance = antigen.newAgent();
    	    instance.setVariable<float>("x", x);
            instance.setVariable<float>("y", y);
            instance.setVariable<float>("vx", 0.0f);
            instance.setVariable<float>("vy", 0.0f);
            instance.setVariable<int>("cell_state", CS_ACTIVE);
    	    instance.setVariable<int>("type", AG_MOLECOLE);
    	    instance.setVariable<float>("mass", 0.1f);
	}
    }
}

int main(int argc, char** argv) {
    // Parse arguments.
    SimulationParams params = parseArguments(argc, argv);

    // Model creation
    flamegpu::ModelDescription model("ImmSimModel");

    /* Initialization of the simulation environment through parameters.
     * The enviroment contains global data, accessible by every agent during
     * the simulation.
     */
    flamegpu::EnvironmentDescription env = model.Environment();
    // Step number.
    env.newProperty<unsigned int>("total_steps", params.steps);
    // Grid dimensions.
    env.newProperty<float>("width", params.width);
    env.newProperty<float>("height", params.height);
    // Cell/Molecule counts.
    env.newProperty<unsigned int>("B_cell_num", params.b_cell_num);
    env.newProperty<unsigned int>("T_cell_num", params.t_cell_num);
    env.newProperty<unsigned int>("antigen_num", params.antigen_num);
    env.newProperty<float>("interaction_radius", PROXIMITY_DIST);

/************************************ AGENT DEFINITION *************************************
 *
 *	Agents are the minimum unit of work inside the FLAMEGPU framewok, each one
 * 	has its own variables.
 */

    // B-Lymphocyte agent.
    flamegpu::AgentDescription bLymph = model.newAgent("B-Lymphocyte");
    bLymph.newVariable<float>("x");
    bLymph.newVariable<float>("y");
    bLymph.newVariable<float>("vx", 0.0f);
    bLymph.newVariable<float>("vy", 0.0f);
    bLymph.newVariable<int>("cell_state", CS_INTERNALIZED);
    bLymph.newVariable<int>("type", B_CELL);
    bLymph.newVariable<float>("mass", 0.2f);
    bLymph.newVariable<uint8_t>("has_interacted", false);
    bLymph.newVariable<unsigned char, RECEPTOR_SIZE>("receptor");

    // T-Lymphocyte agent.
    flamegpu::AgentDescription tLymph = model.newAgent("T-Lymphocyte");
    tLymph.newVariable<float>("x");
    tLymph.newVariable<float>("y");
    tLymph.newVariable<float>("vx", 0.0f);
    tLymph.newVariable<float>("vy", 0.0f);
    tLymph.newVariable<int>("cell_state", CS_ACTIVE);
    tLymph.newVariable<int>("type", T_CELL);
    tLymph.newVariable<float>("mass", 0.2f);
    tLymph.newVariable<uint8_t>("has_interacted", false);

    // Antigen agent.
    flamegpu::AgentDescription antigen = model.newAgent("Antigen");
    antigen.newVariable<float>("x");
    antigen.newVariable<float>("y");
    antigen.newVariable<float>("vx", 0.0f);
    antigen.newVariable<float>("vy", 0.0f);
    antigen.newVariable<int>("cell_state", CS_ACTIVE);
    antigen.newVariable<int>("type", AG_MOLECOLE);
    antigen.newVariable<float>("mass", 0.1f);
    antigen.newVariable<unsigned char, RECEPTOR_SIZE>("epitope");

    // Antibody agent.
    flamegpu::AgentDescription antibody = model.newAgent("Antibody");
    antibody.newVariable<float>("x");
    antibody.newVariable<float>("y");
    antibody.newVariable<float>("vx", 0.0f);
    antibody.newVariable<float>("vy", 0.0f);
    antibody.newVariable<int>("cell_state", CS_ACTIVE);
    antibody.newVariable<int>("type", AB_MOLECOLE);
    antibody.newVariable<uint8_t>("has_interacted", false);
    antibody.newVariable<float>("mass", 0.1f);

/************************************ MESSAGE DEFINITION *************************************
 *
 *	Messages are the medium to exchange information between agents, there are
 *	more types of message other the "Spatial2D", but this one in particular
 *	is highly suitable for this type of simulation.
 */

    // Antigen position message
    flamegpu::MessageSpatial2D::Description antigen_pos = model.newMessage<flamegpu::MessageSpatial2D>("antigen_pos");
    antigen_pos.setMin(0.0f, 0.0f);
    antigen_pos.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    antigen_pos.setRadius(env.getProperty<float>("interaction_radius"));
    antigen_pos.newVariable<flamegpu::id_t>("id");
    antigen_pos.newVariable<unsigned char, RECEPTOR_SIZE>("receptor");

    // Antigen killing message
    flamegpu::MessageSpatial2D::Description antigen_kill = model.newMessage<flamegpu::MessageSpatial2D>("antigen_kill");
    antigen_kill.setMin(0.0f, 0.0f);
    antigen_kill.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    antigen_kill.setRadius(env.getProperty<float>("interaction_radius"));
    antigen_kill.newVariable<flamegpu::id_t>("id");

    // B-Lymphocyte position message
    flamegpu::MessageSpatial2D::Description bcell_pos = model.newMessage<flamegpu::MessageSpatial2D>("bcell_pos");
    bcell_pos.setMin(0.0f, 0.0f);
    bcell_pos.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    bcell_pos.setRadius(env.getProperty<float>("interaction_radius"));
    bcell_pos.newVariable<flamegpu::id_t>("id");
    bcell_pos.newVariable<int>("cell_state");
    bcell_pos.newVariable<uint8_t>("has_interacted");

    // T-Lymphocyte stimulation message
    flamegpu::MessageSpatial2D::Description tcell_stimulate = model.newMessage<flamegpu::MessageSpatial2D>("tcell_stimulate");
    tcell_stimulate.setMin(0.0f, 0.0f);
    tcell_stimulate.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    tcell_stimulate.setRadius(env.getProperty<float>("interaction_radius"));
    tcell_stimulate.newVariable<flamegpu::id_t>("id");

/********************************* BINDING AGENT FUNCTION TO AGENTS **********************************/

    // Antigen agent functions
    antigen.newFunction("antigen_send_position", antigen_send_pos)
        .setMessageOutput("antigen_pos");
    auto recv_kill = antigen.newFunction("antigen_receive_kill", antigen_receive_kill);
    recv_kill.setMessageInput("antigen_kill");
    // Allows agent death inside this function.
    recv_kill.setAllowAgentDeath(true);
    antigen.newFunction("diffuse_entities", diffuse_entities);
    
    // Antibody agent functions
    /* For input/output agent function is necessary to keep track 
     * of the AgentFunctionDescription returned by the newFunction
     * method, since the reference is lost after the first setMessageInput call.
     */
    auto a_i = antibody.newFunction("antibody_interaction", antibody_interaction);
    a_i.setMessageInput("antigen_pos");
    a_i.setMessageOutput("antigen_kill");
    antibody.newFunction("diffuse_entities", diffuse_entities);

    // B-Lymphocyte agent functions
    auto bind = bLymph.newFunction("b_cell_binding", b_cell_binding);
    bind.setMessageInput("antigen_pos");
    bind.setMessageOutput("antigen_kill");
    bLymph.newFunction("b_cell_send_position", b_cell_send_pos)
        .setMessageOutput("bcell_pos");
    bLymph.newFunction("b_cell_receive_stimulus", b_cell_receive_stimulus)
        .setMessageInput("tcell_stimulate");
    auto spawn = bLymph.newFunction("spawn_antibodies", spawn_antibodies);
    spawn.setAgentOutput("Antibody");
    bLymph.newFunction("diffuse_entities", diffuse_entities);

    // T-Lymphocyte agent functions
    auto stim = tLymph.newFunction("t_cell_send_stimulus", t_cell_send_stimulus);
    stim.setMessageInput("bcell_pos");
    stim.setMessageOutput("tcell_stimulate");
    tLymph.newFunction("diffuse_entities", diffuse_entities);

/************************************ LAYER DEFINITION *************************************
 *
 *	Layers define the order execution order of agent functions throughout
 *	the timestep, functions in the same layer execute in parallel.
 *	Each layer can contain ether a host function or more agent functions.
 */

    // Layer creation
    // 1st Layer
    flamegpu::LayerDescription layer1 = model.newLayer("Layer 1");
    layer1.addAgentFunction(antigen.getFunction("antigen_send_position"));

    // 2nd Layer
    flamegpu::LayerDescription layer2 = model.newLayer("Layer 2");
    layer2.addAgentFunction(antibody.getFunction("antibody_interaction"));

    // 3rd Layer
    flamegpu::LayerDescription layer3 = model.newLayer("Layer 3");
    layer3.addAgentFunction(antigen.getFunction("antigen_receive_kill"));

    // 4th Layer
    flamegpu::LayerDescription layer4 = model.newLayer("Layer 4");
    layer4.addAgentFunction(antigen.getFunction("antigen_send_position"));

    // 5th Layer
    flamegpu::LayerDescription layer5 = model.newLayer("Layer 5");
    layer5.addAgentFunction(bLymph.getFunction("b_cell_binding"));

    // 6th Layer
    flamegpu::LayerDescription layer6 = model.newLayer("Layer 6");
    layer6.addAgentFunction(antigen.getFunction("antigen_receive_kill"));

    // 7th Layer
    flamegpu::LayerDescription layer7 = model.newLayer("Layer 7");
    layer7.addAgentFunction(bLymph.getFunction("b_cell_send_position"));

    // 8th Layer
    flamegpu::LayerDescription layer8 = model.newLayer("Layer 8");
    layer8.addAgentFunction(tLymph.getFunction("t_cell_send_stimulus"));

    // 9th Layer
    flamegpu::LayerDescription layer9 = model.newLayer("Layer 9");
    layer9.addAgentFunction(bLymph.getFunction("b_cell_receive_stimulus"));

    // 11th Layer
    flamegpu::LayerDescription layer11 = model.newLayer("Layer 11");
    layer11.addAgentFunction(bLymph.getFunction("spawn_antibodies"));

    // 13th Layer
    flamegpu::LayerDescription layer13 = model.newLayer("Layer 13");
    layer13.addAgentFunction(bLymph.getFunction("diffuse_entities"));
    layer13.addAgentFunction(tLymph.getFunction("diffuse_entities"));
    layer13.addAgentFunction(antibody.getFunction("diffuse_entities"));
    layer13.addAgentFunction(antigen.getFunction("diffuse_entities"));

    // 14th Layer
    flamegpu::LayerDescription layer14 = model.newLayer("Layer 14");
    layer14.addHostFunction(reinsert_antigens);

#ifdef __PLOT__

    // 15th Layer
    flamegpu::LayerDescription layer15 = model.newLayer("Layer 15");
    layer15.addHostFunction(plot_agents);

#endif

    // Set the initialization function defined above.
    model.addInitFunction(init_agents);

    // Set the number of steps.
    flamegpu::CUDASimulation sim(model);
    sim.SimulationConfig().steps = params.steps;

    // Compute execution time.
    auto start = std::chrono::high_resolution_clock::now();
    sim.simulate();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    printf("Simulation time: %f s\n", diff.count());
    return 0;
}
