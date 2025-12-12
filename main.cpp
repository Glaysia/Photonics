
#include "mylib.hpp"
#include "sweep_runner.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace
{

photonics::SimulationSummary dummy_simulation(const photonics::SweepParameters &params)
{
    photonics::SimulationSummary summary{};
    summary.frequency = 0.25 + 0.10 * params.delta_x1; // arbitrary mapping for demo
    summary.quality_factor = 1.0e4 * (1.0 + params.delta_x2);
    summary.mode_volume = 0.9 + 0.5 * (params.r_edge - (params.hole_radius - 0.01 * params.lattice_constant));
    summary.flux_top = 2.0 + params.delta_x1;
    summary.flux_side = 1.0 + params.delta_x2;
    return summary;
}

void run_sweep_runner_self_check()
{
    constexpr double a = 0.42; // µm
    const double r_edge_default = 0.29 * a - 0.01 * a;

    photonics::SweepRunner runner;
    const auto grid = runner.make_grid({0.15 * a, 0.18 * a}, {0.05 * a}, {r_edge_default, r_edge_default - 0.005}, a);
    assert(grid.size() == 4);

    const auto records = runner.run_sweep(grid, dummy_simulation);
    assert(records.size() == grid.size());

    const std::string csv_path = "test_sweep.csv";
    const std::string json_path = "test_sweep.json";
    runner.write_csv(csv_path, records);
    runner.write_json(json_path, records);

    std::ifstream csv_file(csv_path);
    std::string line;
    std::size_t line_count = 0;
    while (std::getline(csv_file, line))
    {
        ++line_count;
    }
    assert(line_count == records.size() + 1); // header + rows

    std::ifstream json_file(json_path);
    const std::string json_text((std::istreambuf_iterator<char>(json_file)), std::istreambuf_iterator<char>());
    assert(json_text.find("\"delta_x1\": 0.063000") != std::string::npos);
    assert(json_text.find("\"quality_factor\"") != std::string::npos);
}

void run_demo_sweep()
{
    constexpr double a = 0.42; // µm
    const double r_edge_default = 0.29 * a - 0.01 * a;

    photonics::SweepRunner runner;
    const auto grid = runner.make_grid({0.15 * a}, {0.05 * a}, {r_edge_default}, a);
    const auto records = runner.run_sweep(grid, dummy_simulation);

    const std::string csv_path = "sample_sweep.csv";
    const std::string json_path = "sample_sweep.json";
    runner.write_csv(csv_path, records);
    runner.write_json(json_path, records);

    std::cout << "[SweepRunner demo] completed " << records.size() << " sweep point(s)\n";
    for (const auto &record : records)
    {
        std::cout << "  dx1=" << record.params.delta_x1
                  << " dx2=" << record.params.delta_x2
                  << " r_edge=" << record.params.r_edge
                  << " -> Q=" << record.summary.quality_factor
                  << " freq=" << record.summary.frequency
                  << " V=" << record.summary.mode_volume << "\n";
    }
    std::cout << "CSV written to " << csv_path << "\n";
    std::cout << "JSON written to " << json_path << "\n";
}

} // namespace

int main()
{
    run_sweep_runner_self_check();
    run_demo_sweep();
    return 0;
}
