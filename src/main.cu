#include <string>
#include <random>
// Not needed python libraries in matplotlib.
#define WITHOUT_NUMPY
#define WITHOUT_OPENCV
#include "matplotlibcpp.h"
#include "math.cuh"
#include "entity.cuh"
#include "flamegpu/flamegpu.h"
// necessary in order to use plot functions.
namespace plt = matplotlibcpp;

#define DEFAULT_WIDTH 500.0f
#define DEFAULT_HEIGHT 500.0f
#define DEFAULT_B_CELL_NUM 100
#define DEFAULT_T_CELL_NUM 50
#define DEFAULT_ANTIGEN_NUM 200
#define DEFAULT_STEPS 5000

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
        // if message is in range and not out of bounds
        const float x2 = msg.getVariable<float>("x");
        const float y2 = msg.getVariable<float>("y");
        if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
            0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
        )
        {
            // check if antigen is in range for interaction.
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
                    unsigned char ag_receptor[RECEPTOR_SIZE];
                    unsigned char receptor[RECEPTOR_SIZE];
		    for (int i=0; i<RECEPTOR_SIZE; i++) {
			ag_receptor[i] = entity.getVariable<unsigned char, RECEPTOR_SIZE>("receptor", i);
			receptor[i] = FLAMEGPU->getVariable<unsigned char, RECEPTOR_SIZE>("receptor", i);

		    }
                    if (can_entities_bind(receptor, ag_receptor)) {
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
            // if message is in range and not out of bounds
            const float x2 = b_cell.getVariable<float>("x");
            const float y2 = b_cell.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // check if antigen is in range for interaction.
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
                        // update interaction flag
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
            // if message is in range and not out of bounds
            const float x2 = stimulus.getVariable<float>("x");
            const float y2 = stimulus.getVariable<float>("y");
            if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
            )
            {
                // check if antigen is in range for interaction.
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

FLAMEGPU_AGENT_FUNCTION(send_position, flamegpu::MessageNone, flamegpu::MessageSpatial2D) {
    const float x = FLAMEGPU->getVariable<float>("x");
    const float y = FLAMEGPU->getVariable<float>("y");
    FLAMEGPU->message_out.setVariable<float>("x", x);
    FLAMEGPU->message_out.setVariable<float>("y", y);
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(spawn_antibodies, flamegpu::MessageSpatial2D, flamegpu::MessageNone) {
    const CellState state = static_cast<CellState>(FLAMEGPU->getVariable<int>("cell_state"));
    const uint8_t has_interacted = FLAMEGPU->getVariable<uint8_t>("has_interacted");
    const int range = 5; // Define a range for spawning antibodies around the B-cell
    int spawned_count = 0;
    if (state == CS_STIMULATED && !has_interacted) {
        const float RADIUS = FLAMEGPU->message_in.radius();
        const float x = FLAMEGPU->getVariable<float>("x");
        const float y = FLAMEGPU->getVariable<float>("y");
        for (int dx = -range; dx <= range && spawned_count < AB_CREATED_PER_CELL; dx++) {
            for (int dy = -range; dy <= range && spawned_count < AB_CREATED_PER_CELL; dy++) {
                if (dx == 0 && dy == 0) continue; // Skip the B-cell's own position
                for (auto msg : FLAMEGPU->message_in(x, y)) {
                    // if message is in range and not out of bounds
                    const float x2 = msg.getVariable<float>("x");
                    const float y2 = msg.getVariable<float>("y");
                    if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                        0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
                    )
                    {
                        // check if antigen is in range for interaction.
                        float x21 = x2 - x;
                        float y21 = y2 - y;
                        const float separation = sqrt( x21*x21 + y21*y21 );
                        if (separation < RADIUS && separation > 0.0f) {
                            if (x2 == (x + dx) && y2 == (y + dy)) {
                                break; // Position occupied, skip to next
                            }
                        }
                    }
                }
                // If position is unoccupied, spawn antibody
		FLAMEGPU->agent_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x") + dx);
		FLAMEGPU->agent_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y") + dy);
                spawned_count++;
            }
        }
        FLAMEGPU->setVariable<uint8_t>("has_interacted", true);
        FLAMEGPU->setVariable<int>("cell_state", CS_INTERNALIZED);
    }
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(diffuse_entities, flamegpu::MessageSpatial2D, flamegpu::MessageNone) {
    const float width = FLAMEGPU->environment.getProperty<float>("width");
    const float height = FLAMEGPU->environment.getProperty<float>("height");
    const float x = FLAMEGPU->getVariable<float>("x");
    const float y = FLAMEGPU->getVariable<float>("y");
    float vx = FLAMEGPU->getVariable<float>("vx");
    float vy = FLAMEGPU->getVariable<float>("vy");
    const double mass = FLAMEGPU->getVariable<float>("mass");

    // Box-Muller
    double r1 = randdouble();
    double r2 = randdouble();
    double random_x = sqrt(-2 * log(r1)) * cos(2 * PI * r2);
    double random_y = sqrt(-2 * log(r1)) * sin(2 * PI * r2);

    // Langevin equation
    vx += langevin(vx, random_x, mass);
    vy += langevin(vy, random_y, mass);

    // Update position, checking boundaries.
    double new_x = floor((x + vx * TIME_FACTOR) > width ? width : (x + vx * TIME_FACTOR < 0) ? 0 : x + vx * TIME_FACTOR);
    double new_y = floor((y + vy * TIME_FACTOR) > height ? height : (y + vy * TIME_FACTOR < 0) ? 0 : y + vy * TIME_FACTOR);

    uint8_t occupied = false;
    const float RADIUS = FLAMEGPU->message_in.radius();
    // Iterate over all received position messages
    for (auto msg : FLAMEGPU->message_in(x, y)) {
        // if message is in range and not out of bounds
        const float x2 = msg.getVariable<float>("x");
        const float y2 = msg.getVariable<float>("y");
        if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
            0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
        )
        {
            // check if antigen is in range for interaction.
            float x21 = x2 - x;
            float y21 = y2 - y;
            const float separation = sqrt( x21*x21 + y21*y21 );
            if (separation < RADIUS && separation > 0.0f) {
                occupied = x2 == new_x && y2 == new_y;
            }
        }
    }
    if (!occupied) {
        FLAMEGPU->setVariable<float>("x", new_x);
        FLAMEGPU->setVariable<float>("y", new_y);
    } else {
        // Search for an unoccupied position adjacent to the new position
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                if (dx == 0 && dy == 0) continue; // Skip the original position
                float adj_x = new_x + dx;
                float adj_y = new_y + dy;
                // Check bounds
                if (adj_x >= 0 && adj_x < width && adj_y >= 0 && adj_y < height) {
                    uint8_t adj_occupied = false;
                    // Re-iterate over messages to check if adjacent position is occupied.
                    for (auto msg : FLAMEGPU->message_in(x, y)) {
                        // if message is in range and not out of bounds
                        const float x2 = msg.getVariable<float>("x");
                        const float y2 = msg.getVariable<float>("y");
                        if ( 0 <= x2 && x2 < FLAMEGPU->environment.getProperty<float>("width") &&
                            0 <= y2 && y2 < FLAMEGPU->environment.getProperty<float>("height")
                        )
                        {
                            // check if antigen is in range for interaction.
                            float x21 = x2 - x;
                            float y21 = y2 - y;
                            const float separation = sqrt( x21*x21 + y21*y21 );
                            if (separation < RADIUS && separation > 0.0f) {
                                if (x2 == adj_x && y2 == adj_y) {
                                    // Whenever it is, skip to next adjacent position.
                                    adj_occupied = true;
                                    break;
                                }
                            }
                        }
                    }
                    // If not occupied, move agent there.
                    if (!adj_occupied) {
                        FLAMEGPU->setVariable<float>("x", adj_x);
                        FLAMEGPU->setVariable<float>("y", adj_y);
                        // If a valid position is found, exit both loops.
                        occupied = false;
                        break;
                    }
                }
            }
            if (!occupied) break;
        }
    }
    // Update velocity variables.
    FLAMEGPU->setVariable<float>("vx", vx);
    FLAMEGPU->setVariable<float>("vy", vy);
    // Reset interaction flag.
    EntityType type = static_cast<EntityType>(FLAMEGPU->getVariable<int>("type"));
    if (type != AG_MOLECOLE){
        FLAMEGPU->setVariable<uint8_t>("has_interacted", false);
    }
    // If no valid position is found, the agent remains in its current position.
    return flamegpu::ALIVE;
}

FLAMEGPU_INIT_FUNCTION(init_agents) {
    const unsigned int width = FLAMEGPU->environment.getProperty<float>("width");
    const unsigned int height = FLAMEGPU->environment.getProperty<float>("height");
    const unsigned int b_cell_num = FLAMEGPU->environment.getProperty<unsigned int>("B_cell_num");
    const unsigned int t_cell_num = FLAMEGPU->environment.getProperty<unsigned int>("T_cell_num");
    const unsigned int antigen_num = FLAMEGPU->environment.getProperty<unsigned int>("antigen_num");
    const unsigned int total_agents = b_cell_num + t_cell_num + antigen_num;
    // create rng
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, MAX_ENTITYTYPE - 1); // define the range, excluding antibodies

    auto b_cell = FLAMEGPU->agent("B-Lymphocyte");
    auto t_cell = FLAMEGPU->agent("T-Lymphocyte");
    auto antigen = FLAMEGPU->agent("Antigen");
    auto antibody = FLAMEGPU->agent("Antibody");

    unsigned int counter = 0, b = 0, t = 0, ag = 0;

    // How many agents in the grid
    unsigned int cols = static_cast<unsigned int>(std::ceil(std::sqrt(total_agents)));
    unsigned int rows = static_cast<unsigned int>(std::ceil(static_cast<float>(total_agents) / cols));

    float cell_width = static_cast<float>(width) / cols;
    float cell_height = static_cast<float>(height) / rows;

    for (unsigned int i = 0; i < rows && counter < total_agents; ++i) {
        for (unsigned int j = 0; j < cols && counter < total_agents; ++j, counter++) {
            unsigned int type;
            uint8_t full = true;
            while(full && !(ag >= antigen_num || b >= b_cell_num || t >= t_cell_num)) {
                // Choose a random type.
                type = distr(gen);
                switch (static_cast<EntityType>(type)) {
                    case B_CELL:
                        full = b >= b_cell_num;
                        break;
                    case T_CELL:
                        full = t >= t_cell_num;
                        break;
                    case AG_MOLECOLE:
                        full = ag >= antigen_num;
                        break;
                    default:
                        break;
                }
            }
            switch(static_cast<EntityType>(type)) {
                case B_CELL: {
		    auto instance = b_cell.newAgent();
		    b++;
                    instance.setVariable<int>("cell_state", CS_INTERNALIZED);
                    instance.setVariable<int>("type", B_CELL);
                    instance.setVariable<float>("mass", 0.2f);
                    instance.setVariable<uint8_t>("has_interacted", false);
                    for (unsigned int k = 0; k < RECEPTOR_SIZE; ++k) {
                        instance.setVariable<unsigned char, RECEPTOR_SIZE>("receptor", k, randbyte());
                    }
		    instance.setVariable<float>("x", j);
            	    instance.setVariable<float>("y", i);
            	    instance.setVariable<float>("vx", 0.0f);
            	    instance.setVariable<float>("vy", 0.0f);
                    break;
                }
                case T_CELL: {
		    auto instance = t_cell.newAgent();
		    t++;
                    instance.setVariable<int>("cell_state", CS_ACTIVE);
                    instance.setVariable<int>("type", T_CELL);
                    instance.setVariable<float>("mass", 0.2f);
                    instance.setVariable<uint8_t>("has_interacted", false);
		    instance.setVariable<float>("x", j);
            	    instance.setVariable<float>("y", i);
            	    instance.setVariable<float>("vx", 0.0f);
            	    instance.setVariable<float>("vy", 0.0f);
                    break;
                }
                case AG_MOLECOLE: {
		    auto instance = antigen.newAgent();
		    ag++;
                    instance.setVariable<int>("cell_state", CS_ACTIVE);
                    instance.setVariable<int>("type", AG_MOLECOLE);
                    instance.setVariable<float>("mass", 0.1f);
                    for (unsigned int k = 0; k < RECEPTOR_SIZE; ++k) {
                        instance.setVariable<unsigned char, RECEPTOR_SIZE>("epitope", k, randbyte());
                    }
		    instance.setVariable<float>("x", j);
            	    instance.setVariable<float>("y", i);
            	    instance.setVariable<float>("vx", 0.0f);
            	    instance.setVariable<float>("vy", 0.0f);
                    break;
                }
                default:
                    break;
            }
        }
    }
}

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
        auto b_cells = std::move(FLAMEGPU->agent("B-Lymphocyte").getPopulationData());
        auto t_cells = std::move(FLAMEGPU->agent("T-Lymphocyte").getPopulationData());
        auto antigens = std::move(FLAMEGPU->agent("Antigen").getPopulationData());
        auto antibodies = std::move(FLAMEGPU->agent("Antibody").getPopulationData());

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

int main(int argc, char** argv) {
    float width, height;
    unsigned int b_cell_num, t_cell_num, antigen_num, steps;
    if (argc == 7) {
        width = std::stof(argv[1]);
        height = std::stof(argv[2]);
        b_cell_num = std::stoi(argv[3]);
        t_cell_num = std::stoi(argv[4]);
        antigen_num = std::stoi(argv[5]);
        steps = std::stoi(argv[6]);
    } else {
        width = DEFAULT_WIDTH;
        height = DEFAULT_HEIGHT;
        b_cell_num = DEFAULT_B_CELL_NUM;
        t_cell_num = DEFAULT_T_CELL_NUM;
        antigen_num = DEFAULT_ANTIGEN_NUM;
        steps = DEFAULT_STEPS;
    }

    // Model creation
    flamegpu::ModelDescription model("ImmSimModel");

    // Environment.
    flamegpu::EnvironmentDescription env = model.Environment();
    // Grid dimensions.
    env.newProperty<float>("width", width);
    env.newProperty<float>("height", height);
    // Cell/Molecule counts.
    env.newProperty<unsigned int>("B_cell_num", b_cell_num);
    env.newProperty<unsigned int>("T_cell_num", t_cell_num);
    env.newProperty<unsigned int>("antigen_num", antigen_num);
    // Constant radius to determine if two agents can interact.
    env.newProperty<float>("interaction_radius", PROXIMITY_DIST);
    // Constant radius for spatial messages.
    env.newProperty<float>("spatial_radius", 10.0f);

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
    // antibody.newState("default");

    // Antigen position message
    flamegpu::MessageSpatial2D::Description antigen_pos = model.newMessage<flamegpu::MessageSpatial2D>("antigen_pos");
    antigen_pos.setMin(0.0f, 0.0f);
    antigen_pos.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    antigen_pos.setRadius(env.getProperty<float>("interaction_radius"));
    //antigen_pos.newVariable<float>("x");
    //antigen_pos.newVariable<float>("y");
    antigen_pos.newVariable<flamegpu::id_t>("id");
    antigen_pos.newVariable<unsigned char, RECEPTOR_SIZE>("receptor");

    // Antigen killing message
    flamegpu::MessageSpatial2D::Description antigen_kill = model.newMessage<flamegpu::MessageSpatial2D>("antigen_kill");
    antigen_kill.setMin(0.0f, 0.0f);
    antigen_kill.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    antigen_kill.setRadius(env.getProperty<float>("interaction_radius"));
    //antigen_kill.newVariable<float>("x");
    //antigen_kill.newVariable<float>("y");
    antigen_kill.newVariable<flamegpu::id_t>("id");

    // B-Lymphocyte position message
    flamegpu::MessageSpatial2D::Description bcell_pos = model.newMessage<flamegpu::MessageSpatial2D>("bcell_pos");
    bcell_pos.setMin(0.0f, 0.0f);
    bcell_pos.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    bcell_pos.setRadius(env.getProperty<float>("interaction_radius"));
    //bcell_pos.newVariable<float>("x");
    //bcell_pos.newVariable<float>("y");
    bcell_pos.newVariable<flamegpu::id_t>("id");
    bcell_pos.newVariable<int>("cell_state");
    bcell_pos.newVariable<uint8_t>("has_interacted");

    // T-Lymphocyte stimulation message
    flamegpu::MessageSpatial2D::Description tcell_stimulate = model.newMessage<flamegpu::MessageSpatial2D>("tcell_stimulate");
    tcell_stimulate.setMin(0.0f, 0.0f);
    tcell_stimulate.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    tcell_stimulate.setRadius(env.getProperty<float>("interaction_radius"));
    //tcell_stimulate.newVariable<float>("x");
    //tcell_stimulate.newVariable<float>("y");
    tcell_stimulate.newVariable<flamegpu::id_t>("id");

    // Generic agent position message for diffusion
    flamegpu::MessageSpatial2D::Description position = model.newMessage<flamegpu::MessageSpatial2D>("position");
    position.setMin(0.0f, 0.0f);
    position.setMax(env.getProperty<float>("width"), env.getProperty<float>("height"));
    position.setRadius(env.getProperty<float>("spatial_radius"));
    //position.newVariable<float>("x");
    //position.newVariable<float>("y");

    // Antigen agent functions
    antigen.newFunction("antigen_send_position", antigen_send_pos)
        .setMessageOutput("antigen_pos");
    auto recv_kill = antigen.newFunction("antigen_receive_kill", antigen_receive_kill);
    recv_kill.setMessageInput("antigen_kill");
    // Allows agent death inside this function.
    recv_kill.setAllowAgentDeath(true);
    antigen.newFunction("send_position", send_position)
        .setMessageOutput("position");
    antigen.newFunction("diffuse_entities", diffuse_entities)
        .setMessageInput("position");
    
    // Antibody agent functions
    /* For input/output agent function is necessary to keep track 
     * of the AgentFunctionDescription returned by the newFunction
     * method, since the reference is lost after the first setMessageInput call.
     */
    auto a_i = antibody.newFunction("antibody_interaction", antibody_interaction);
    a_i.setMessageInput("antigen_pos");
    a_i.setMessageOutput("antigen_kill");
    antibody.newFunction("send_position", send_position)
        .setMessageOutput("position");
    antibody.newFunction("diffuse_entities", diffuse_entities)
        .setMessageInput("position");

    // B-Lymphocyte agent functions
    auto bind = bLymph.newFunction("b_cell_binding", b_cell_binding);
    bind.setMessageInput("antigen_pos");
    bind.setMessageOutput("antigen_kill");
    bLymph.newFunction("b_cell_send_position", b_cell_send_pos)
        .setMessageOutput("bcell_pos");
    bLymph.newFunction("b_cell_receive_stimulus", b_cell_receive_stimulus)
        .setMessageInput("tcell_stimulate");
    auto spawn = bLymph.newFunction("spawn_antibodies", spawn_antibodies);
    spawn.setMessageInput("position");
    spawn.setAgentOutput("Antibody");
    bLymph.newFunction("send_position", send_position)
        .setMessageOutput("position");
    bLymph.newFunction("diffuse_entities", diffuse_entities)
        .setMessageInput("position");

    // T-Lymphocyte agent functions
    auto stim = tLymph.newFunction("t_cell_send_stimulus", t_cell_send_stimulus);
    stim.setMessageInput("bcell_pos");
    stim.setMessageOutput("tcell_stimulate");
    tLymph.newFunction("send_position", send_position)
        .setMessageOutput("position");
    tLymph.newFunction("diffuse_entities", diffuse_entities)
        .setMessageInput("position");

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

    /* 10th Layer
    flamegpu::LayerDescription layer10 = model.newLayer("Layer 10");
    layer10.addAgentFunction(bLymph.getFunction("send_position"));
    layer10.addAgentFunction(tLymph.getFunction("send_position"));
    layer10.addAgentFunction(antibody.getFunction("send_position"));
    layer10.addAgentFunction(antigen.getFunction("send_position"));
    */
    auto layer10_b = model.newLayer("Layer 10 B-Lymphocyte");
    layer10_b.addAgentFunction(bLymph.getFunction("send_position"));

    auto layer10_t = model.newLayer("Layer 10 T-Lymphocyte");
    layer10_t.addAgentFunction(tLymph.getFunction("send_position"));

    auto layer10_a = model.newLayer("Layer 10 Antibody");
    layer10_a.addAgentFunction(antibody.getFunction("send_position"));

    auto layer10_ag = model.newLayer("Layer 10 Antigen");
    layer10_ag.addAgentFunction(antigen.getFunction("send_position"));


    // 11th Layer
    flamegpu::LayerDescription layer11 = model.newLayer("Layer 11");
    layer11.addAgentFunction(bLymph.getFunction("spawn_antibodies"));

    /* 12th Layer
    flamegpu::LayerDescription layer12 = model.newLayer("Layer 12");
    layer12.addAgentFunction(bLymph.getFunction("send_position"));
    layer12.addAgentFunction(tLymph.getFunction("send_position"));
    layer12.addAgentFunction(antibody.getFunction("send_position"));
    layer12.addAgentFunction(antigen.getFunction("send_position"));
    */
    auto layer12_b = model.newLayer("Layer 12 B-Lymphocyte");
    layer12_b.addAgentFunction(bLymph.getFunction("send_position"));

    auto layer12_t = model.newLayer("Layer 12 T-Lymphocyte");
    layer12_t.addAgentFunction(tLymph.getFunction("send_position"));

    auto layer12_a = model.newLayer("Layer 12 Antibody");
    layer12_a.addAgentFunction(antibody.getFunction("send_position"));

    auto layer12_ag = model.newLayer("Layer 12 Antigen");
    layer12_ag.addAgentFunction(antigen.getFunction("send_position"));

    // 13th Layer
    flamegpu::LayerDescription layer13 = model.newLayer("Layer 13");
    layer13.addAgentFunction(bLymph.getFunction("diffuse_entities"));
    layer13.addAgentFunction(tLymph.getFunction("diffuse_entities"));
    layer13.addAgentFunction(antibody.getFunction("diffuse_entities"));
    layer13.addAgentFunction(antigen.getFunction("diffuse_entities"));

    // 14th Layer
    flamegpu::LayerDescription layer14 = model.newLayer("Layer 14");
    layer14.addHostFunction(plot_agents);

    // Initialization function
    model.addInitFunction(init_agents);

    flamegpu::CUDASimulation sim(model);
    sim.SimulationConfig().steps = steps;
    sim.simulate();
    return 0;
}
