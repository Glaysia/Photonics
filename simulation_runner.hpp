#ifndef PHOTONICS_SIMULATION_RUNNER_HPP
#define PHOTONICS_SIMULATION_RUNNER_HPP

#include "simulation_config.hpp"
#include "source_manager.hpp"

#include <meep.hpp>

#include <complex>
#include <functional>
#include <memory>
#include <vector>

namespace photonics
{

struct TimeSeriesConfig
{
    bool enabled = false;
    meep::vec point{0.0, 0.0, 0.0};
    meep::component component = meep::Ey;
    int sample_every_steps = 10;          // must be >= 1 for uniform sampling
    double start_after_sources = 0.0;     // time in units of a/c after sources end
};

struct SimulationRunConfig
{
    double until_after_sources = 200.0; // time in units of a/c after sources end
    double progress_interval_seconds = 7.0; // wall-clock interval for progress reporting; <=0 disables
    TimeSeriesConfig time_series{};
};

struct SimulationRunResult
{
    double last_source_time = 0.0;
    double stop_time = 0.0;
    int steps = 0;

    double timestep = 0.0;
    double sample_dt = 0.0;
    std::vector<std::complex<double>> samples;

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
