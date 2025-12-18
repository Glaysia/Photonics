
#include "lattice_geometry.hpp"
#include "geometry_export.hpp"
#include "monitor_suite.hpp"
#include "q_analyzer.hpp"
#include "simulation_config.hpp"
#include "source_manager.hpp"
#include "sweep_runner.hpp"

#include <meep.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

namespace
{
void banner(const std::string &title)
{
    std::cout << "\n== " << title << " ==\n";
}

std::string tuple_to_string(const std::tuple<int, int, int> &t)
{
    std::ostringstream oss;
    oss << std::get<0>(t) << " × " << std::get<1>(t) << " × " << std::get<2>(t);
    return oss.str();
}

bool has_flag(int argc, char **argv, const std::string &flag)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == flag)
        {
            return true;
        }
    }
    return false;
}

std::string arg_value(int argc, char **argv, const std::string &flag)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == flag)
        {
            if (i + 1 < argc && argv[i + 1][0] != '-')
            {
                return argv[i + 1];
            }
            return {};
        }
    }
    return {};
}

std::filesystem::path ensure_parent_dir(std::filesystem::path path)
{
    if (path.has_parent_path())
    {
        std::filesystem::create_directories(path.parent_path());
    }
    return path;
}
} // namespace

int main(int argc, char **argv)
{
    meep::initialize mpi(argc, argv);

    // Simulation parameters derived from the Akahane L3 cavity defaults.
    photonics::SimulationParameters sim_params{};
    sim_params.lattice_constant_um = 0.42;
    sim_params.resolution_px_per_a = 20.0;
    sim_params.cell_size_x_in_a = 10.0;
    sim_params.cell_size_y_in_a = 8.0;
    sim_params.cell_size_z_in_a = 5.0;
    sim_params.pml_thickness_in_a = 1.0;
    sim_params.courant = 0.5;

    photonics::SimulationConfig sim_config(sim_params);
    const auto grid_points = sim_config.expected_grid_points();

    banner("Grid / Volume");
    std::cout << "Cell (a units): " << sim_config.cell_size_x_in_a() << " x " << sim_config.cell_size_y_in_a()
              << " x " << sim_config.cell_size_z_in_a() << "\n";
    std::cout << "Expected grid points: " << tuple_to_string(grid_points) << "\n";
    std::cout << "Resolution: " << sim_config.resolution_px_per_a() << " px/a\n";
    std::cout << "PML thickness: " << sim_config.pml_thickness_in_a() << " a\n";
    std::cout << "Courant: " << sim_config.courant() << " (dt=" << sim_config.timestep() << ")\n";

    // Lattice + geometry for L3 defect.
    photonics::LatticeGeometryParams geom_params{};
    geom_params.lattice_constant = sim_params.lattice_constant_um;
    geom_params.hole_radius = 0.29 * geom_params.lattice_constant;
    geom_params.slab_thickness = 0.6 * geom_params.lattice_constant;
    geom_params.delta_x1 = 0.15 * geom_params.lattice_constant;
    geom_params.delta_x2 = 0.05 * geom_params.lattice_constant;
    geom_params.edge_radius = geom_params.hole_radius - 0.01 * geom_params.lattice_constant;

    photonics::LatticeGeometry geometry_builder(geom_params);
    const auto holes = geometry_builder.generate_holes();
    const auto geometry = geometry_builder.build_geometry();

    banner("Geometry");
    std::cout << "Holes (after defect removal): " << holes.size() << "\n";
    std::cout << "Sample hole: center=(" << holes.front().center.x() << ", " << holes.front().center.y()
              << "), r=" << holes.front().radius << "\n";

    if (has_flag(argc, argv, "--export-geometry"))
    {
        const std::string basename = [&]() {
            const auto v = arg_value(argc, argv, "--export-geometry");
            return v.empty() ? std::string("out/l3_geometry") : v;
        }();

        photonics::GeometryMetadata meta{};
        meta.params = geom_params;
        meta.holes = holes;
        meta.units = "um";

        const auto csv_path = ensure_parent_dir(basename + ".csv");
        const auto json_path = ensure_parent_dir(basename + ".json");
        const auto svg_path = ensure_parent_dir(basename + ".svg");

        photonics::write_holes_csv(csv_path.string(), holes);
        photonics::write_geometry_json(json_path.string(), meta);
        photonics::write_holes_svg(svg_path.string(), holes);

        std::cout << "Exported geometry:\n";
        std::cout << "  - " << csv_path << "\n";
        std::cout << "  - " << json_path << "\n";
        std::cout << "  - " << svg_path << "\n";

        if (has_flag(argc, argv, "--export-only"))
        {
            return EXIT_SUCCESS;
        }
    }

    // Structure + sources.
    meep::structure structure = sim_config.make_structure(geometry);
    photonics::SourceManager source_manager;
    const auto sources = source_manager.build_sources();

    banner("Sources");
    std::cout << "Configured sources: " << sources.size() << " (components × positions)\n";
    std::cout << "Bandwidth: " << source_manager.config().bandwidth << " (1/a)\n";

    // Monitoring scaffold.
    photonics::MonitorConfig monitor_cfg{};
    monitor_cfg.cell = {sim_params.cell_size_x_in_a, sim_params.cell_size_y_in_a, sim_params.cell_size_z_in_a, meep::D3};
    monitor_cfg.pml_thickness = sim_params.pml_thickness_in_a;
    monitor_cfg.monitor_thickness = 0.05;

    photonics::MonitorSuite monitors(monitor_cfg);
    monitors.setup_default_flux_boxes();
    monitors.setup_default_field_snapshot(meep::Z, 0.0);
    monitors.add_harminv(meep::vec(0.0, 0.0, 0.0), 0.26, 0.08, meep::Ez);

    banner("Monitors");
    std::cout << "Flux monitors: " << monitors.flux_monitors().size() << "\n";
    std::cout << "Field snapshots: " << monitors.field_monitors().size() << "\n";
    std::cout << "Harminv probes: " << monitors.harminv_monitors().size() << "\n";

    // Dummy Harminv result to demonstrate Q extraction pipeline.
    photonics::QAnalyzer analyzer;
    const std::vector<photonics::HarminvMode> dummy_modes = {
        {0.26, 2.0e-4, 1.0},
        {0.24, 5.0e-4, 0.2},
    };
    const auto qres = analyzer.compute_quality(dummy_modes);
    const auto flux = analyzer.flux_components(1.0, 0.8, 0.1);

    banner("Q Analysis (dummy data)");
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Resonance freq: " << qres.resonance_frequency << "  decay: " << qres.decay_rate
              << "  Q: " << qres.quality_factor << "\n";
    std::cout << "Flux breakdown (top/bottom/lateral): " << flux.upward << " / " << flux.downward << " / "
              << flux.lateral << "\n";

    // Parameter sweep scaffold.
    photonics::SweepRunner sweeper;
    const auto grid = sweeper.make_grid({geom_params.delta_x1}, {geom_params.delta_x2}, {geom_params.edge_radius},
                                        geom_params.lattice_constant);

    banner("Sweep Grid");
    std::cout << "Sweep combinations: " << grid.size() << " (dx1/dx2/r_edge)\n";
    for (const auto &p : grid)
    {
        std::cout << "  dx1=" << p.delta_x1 << ", dx2=" << p.delta_x2 << ", r_edge=" << p.r_edge << "\n";
    }

    banner("Next Steps");
    std::cout << "- Replace dummy Harminv data with actual meep::harminv analysis on fields.\n";
    std::cout << "- Add flux collection using MonitorSuite definitions (meep flux boxes) and feed into QAnalyzer.\n";
    std::cout << "- Drive SweepRunner with a simulation callback that returns SimulationSummary per parameter set.\n";

    return EXIT_SUCCESS;
}
