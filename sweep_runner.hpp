// SweepRunner orchestrates parameter sweeps and result persistence.
#ifndef PHOTONICS_SWEEP_RUNNER_HPP
#define PHOTONICS_SWEEP_RUNNER_HPP

#include <functional>
#include <string>
#include <vector>

namespace photonics
{

struct SweepParameters
{
    double lattice_constant = 1.0;          // a (simulation units)
    double delta_x1 = 0.15 * lattice_constant;
    double delta_x2 = 0.05 * lattice_constant;
    double hole_radius = 0.29 * lattice_constant;
    double r_edge = hole_radius - 0.01 * lattice_constant;
};

struct SimulationSummary
{
    double frequency = 0.0;
    double quality_factor = 0.0;
    double mode_volume = 0.0;
    double flux_top = 0.0;
    double flux_side = 0.0;
};

struct SweepRecord
{
    SweepParameters params{};
    SimulationSummary summary{};
};

class SweepRunner
{
public:
    using SimulationCallback = std::function<SimulationSummary(const SweepParameters &)>;

    SweepRunner() = default;

    std::vector<SweepParameters> make_grid(const std::vector<double> &delta_x1_values,
                                           const std::vector<double> &delta_x2_values,
                                           const std::vector<double> &r_edge_values,
                                           double lattice_constant = 0.42) const;

    std::vector<SweepRecord> run_sweep(const std::vector<SweepParameters> &grid,
                                       const SimulationCallback &callback) const;

    void write_csv(const std::string &path, const std::vector<SweepRecord> &records) const;
    void write_json(const std::string &path, const std::vector<SweepRecord> &records) const;
};

} // namespace photonics

#endif // PHOTONICS_SWEEP_RUNNER_HPP
