#include "sweep_runner.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace photonics
{

std::vector<SweepParameters> SweepRunner::make_grid(const std::vector<double> &delta_x1_values,
                                                    const std::vector<double> &delta_x2_values,
                                                    const std::vector<double> &r_edge_values,
                                                    double lattice_constant) const
{
    std::vector<SweepParameters> grid;
    grid.reserve(delta_x1_values.size() * delta_x2_values.size() * r_edge_values.size());

    for (const auto dx1 : delta_x1_values)
    {
        for (const auto dx2 : delta_x2_values)
        {
            for (const auto r_edge : r_edge_values)
            {
                SweepParameters params;
                params.lattice_constant = lattice_constant;
                params.hole_radius = 0.29 * lattice_constant;
                params.delta_x1 = dx1;
                params.delta_x2 = dx2;
                params.r_edge = r_edge;
                grid.push_back(params);
            }
        }
    }

    return grid;
}

std::vector<SweepRecord> SweepRunner::run_sweep(const std::vector<SweepParameters> &grid,
                                                const SimulationCallback &callback) const
{
    std::vector<SweepRecord> results;
    results.reserve(grid.size());

    for (const auto &params : grid)
    {
        SweepRecord record;
        record.params = params;
        record.summary = callback(params);
        results.push_back(record);
    }

    return results;
}

void SweepRunner::write_csv(const std::string &path, const std::vector<SweepRecord> &records) const
{
    std::ofstream out(path, std::ios::trunc);
    out << "delta_x1,delta_x2,r_edge,frequency,quality_factor,mode_volume,flux_top,flux_side\n";
    out << std::fixed << std::setprecision(6);
    for (const auto &record : records)
    {
        out << record.params.delta_x1 << ","
            << record.params.delta_x2 << ","
            << record.params.r_edge << ","
            << record.summary.frequency << ","
            << record.summary.quality_factor << ","
            << record.summary.mode_volume << ","
            << record.summary.flux_top << ","
            << record.summary.flux_side << "\n";
    }
}

void SweepRunner::write_json(const std::string &path, const std::vector<SweepRecord> &records) const
{
    std::ofstream out(path, std::ios::trunc);
    out << std::fixed << std::setprecision(6);
    out << "[\n";

    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const auto &record = records[i];
        out << "  {\n";
        out << "    \"delta_x1\": " << record.params.delta_x1 << ",\n";
        out << "    \"delta_x2\": " << record.params.delta_x2 << ",\n";
        out << "    \"r_edge\": " << record.params.r_edge << ",\n";
        out << "    \"frequency\": " << record.summary.frequency << ",\n";
        out << "    \"quality_factor\": " << record.summary.quality_factor << ",\n";
        out << "    \"mode_volume\": " << record.summary.mode_volume << ",\n";
        out << "    \"flux_top\": " << record.summary.flux_top << ",\n";
        out << "    \"flux_side\": " << record.summary.flux_side << "\n";
        out << "  }";
        if (i + 1 != records.size())
        {
            out << ",";
        }
        out << "\n";
    }

    out << "]\n";
}

} // namespace photonics
