#ifndef PHOTONICS_SIMULATION_RUNNER_HPP
#define PHOTONICS_SIMULATION_RUNNER_HPP

#include "simulation_config.hpp"
#include "source_manager.hpp"

#include <meep.hpp>

#include <memory>
#include <vector>

namespace photonics
{

struct SimulationRunConfig
{
    double until_after_sources = 200.0; // time in units of a/c after sources end
};

struct SimulationRunResult
{
    double last_source_time = 0.0;
    double stop_time = 0.0;
    int steps = 0;

    std::unique_ptr<meep::structure> structure;
    std::unique_ptr<meep::fields> fields;
};

class SimulationRunner
{
public:
    SimulationRunner() = default;

    SimulationRunResult run(const SimulationConfig &config,
                            const std::vector<meep::geometric_object> &geometry,
                            const std::vector<meep::source> &sources,
                            const SimulationRunConfig &run_config = {}) const;
};

} // namespace photonics

#endif // PHOTONICS_SIMULATION_RUNNER_HPP

