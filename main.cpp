
#include "lattice_geometry.hpp"
#include "geometry_export.hpp"
#include "monitor_suite.hpp"
#include "q_analyzer.hpp"
#include "harminv_runner.hpp"
#include "simulation_config.hpp"
#include "simulation_runner.hpp"
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
            const auto is_number = [](const char *s) {
                if (!s || *s == '\0')
                {
                    return false;
                }
                char *end = nullptr;
                std::strtod(s, &end);
                return end && *end == '\0';
            };

            if (i + 1 < argc && (argv[i + 1][0] != '-' || is_number(argv[i + 1])))
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

int parse_int_env(const char *name)
{
    const char *v = std::getenv(name);
    if (!v || *v == '\0')
    {
        return -1;
    }
    char *end = nullptr;
    const long parsed = std::strtol(v, &end, 10);
    if (!end || *end != '\0')
    {
        return -1;
    }
    return static_cast<int>(parsed);
}

int process_rank()
{
    // Prefer Meep MPI rank if available, otherwise fall back to common MPI env vars.
    if (meep::with_mpi())
    {
        return meep::my_rank();
    }
    for (const char *name : {"OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_PROCID", "MV2_COMM_WORLD_RANK"})
    {
        const int v = parse_int_env(name);
        if (v >= 0)
        {
            return v;
        }
    }
    return 0;
}

void discard_meep_output(const char *)
{
}
} // namespace

int main(int argc, char **argv)
{
    meep::initialize mpi(argc, argv);
    const bool is_master = (process_rank() == 0);
    if (!is_master)
    {
        meep::verbosity = 0;
        meep::set_meep_printf_callback(&discard_meep_output);
        meep::set_meep_printf_stderr_callback(&discard_meep_output);
    }

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

    if (is_master)
    {
        banner("Grid / Volume");
        std::cout << "Cell (a units): " << sim_config.cell_size_x_in_a() << " x " << sim_config.cell_size_y_in_a()
                  << " x " << sim_config.cell_size_z_in_a() << "\n";
        std::cout << "Expected grid points: " << tuple_to_string(grid_points) << "\n";
        std::cout << "Resolution: " << sim_config.resolution_px_per_a() << " px/a\n";
        std::cout << "PML thickness: " << sim_config.pml_thickness_in_a() << " a\n";
        std::cout << "Courant: " << sim_config.courant() << " (dt=" << sim_config.timestep() << ")\n";
    }

    // Lattice + geometry for L3 defect.
    photonics::LatticeGeometryParams geom_params{};
    // Use normalized simulation units (a = 1). Physical scaling is applied after the fact.
    geom_params.lattice_constant = 1.0;
    geom_params.hole_radius = 0.29 * geom_params.lattice_constant;
    geom_params.slab_thickness = 0.6 * geom_params.lattice_constant;
    geom_params.delta_x1 = 0.15 * geom_params.lattice_constant;
    geom_params.delta_x2 = 0.05 * geom_params.lattice_constant;
    geom_params.edge_radius = geom_params.hole_radius - 0.01 * geom_params.lattice_constant;

    photonics::LatticeGeometry geometry_builder(geom_params);
    const auto holes = geometry_builder.generate_holes();
    const auto geometry = geometry_builder.build_geometry();

    if (is_master)
    {
        banner("Geometry");
        std::cout << "Holes (after defect removal): " << holes.size() << "\n";
        std::cout << "Sample hole: center=(" << holes.front().center.x() << ", " << holes.front().center.y()
                  << "), r=" << holes.front().radius << "\n";
    }

    if (is_master && has_flag(argc, argv, "--export-geometry"))
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

    // Sources.
    photonics::SourceManager source_manager;
    const auto sources = source_manager.build_sources();

    if (is_master)
    {
        banner("Sources");
        std::cout << "Configured sources: " << sources.size() << " (components × positions)\n";
        std::cout << "Bandwidth: " << source_manager.config().bandwidth << " (1/a)\n";
    }

    const bool measure_q = has_flag(argc, argv, "--measure-q");
    if (has_flag(argc, argv, "--run-until-after-sources") || measure_q)
    {
        photonics::SimulationRunner runner;
        photonics::SimulationRunConfig run_cfg{};
        run_cfg.until_after_sources = measure_q ? 50000.0 : 200.0;
        if (const auto v = arg_value(argc, argv, "--run-until-after-sources"); !v.empty())
        {
            run_cfg.until_after_sources = std::stod(v);
        }

        if (measure_q)
        {
            run_cfg.time_series.enabled = true;
            run_cfg.time_series.point = meep::vec(0.0, 0.0, 0.0);
            run_cfg.time_series.component = meep::Ey; // TE-like cavity mode is typically in-plane E
            run_cfg.time_series.sample_every_steps = 10;
            run_cfg.time_series.start_after_sources = 0.0;
        }

        const auto res = runner.run(sim_config, geometry, sources, run_cfg);
        if (is_master)
        {
            banner("Run");
            std::cout << "Last source time: " << res.last_source_time << "\n";
            std::cout << "Stop time: " << res.stop_time << "\n";
            std::cout << "Steps: " << res.steps << "\n";
            if (measure_q)
            {
                std::cout << "Sample dt: " << res.sample_dt << " (time units a/c)\n";
                std::cout << "Samples: " << res.samples.size() << "\n";
            }
        }

        if (measure_q)
        {
            photonics::HarminvRunConfig hcfg{};
            const double fc = source_manager.config().center_frequency;
            const double bw = std::max(source_manager.config().bandwidth, 1e-6);
            hcfg.fmin = fc - 2.0 * bw;
            hcfg.fmax = fc + 2.0 * bw;
            hcfg.maxbands = 40;
            hcfg.q_threshold = 100.0;

            const auto modes = photonics::run_harminv(res.samples, res.sample_dt, hcfg);
            photonics::QAnalyzer analyzer;
            const auto qres = analyzer.compute_quality(modes);

            if (is_master)
            {
                banner("Harminv / Q");
                std::cout << std::fixed << std::setprecision(9);
                std::cout << "Search band: [" << hcfg.fmin << ", " << hcfg.fmax << "] (1/a)\n";
                std::cout << "Modes found: " << modes.size() << "\n";
                std::cout << "Resonance freq: " << qres.resonance_frequency << "  decay: " << qres.decay_rate
                          << "  Q: " << qres.quality_factor << "\n";
            }
        }
    }

    // Monitoring scaffold.
    photonics::MonitorConfig monitor_cfg{};
    monitor_cfg.cell = {sim_params.cell_size_x_in_a, sim_params.cell_size_y_in_a, sim_params.cell_size_z_in_a, meep::D3};
    monitor_cfg.pml_thickness = sim_params.pml_thickness_in_a;
    monitor_cfg.monitor_thickness = 0.05;

    photonics::MonitorSuite monitors(monitor_cfg);
    monitors.setup_default_flux_boxes();
    monitors.setup_default_field_snapshot(meep::Z, 0.0);
    monitors.add_harminv(meep::vec(0.0, 0.0, 0.0), 0.26, 0.08, meep::Ez);

    if (is_master)
    {
        banner("Monitors");
        std::cout << "Flux monitors: " << monitors.flux_monitors().size() << "\n";
        std::cout << "Field snapshots: " << monitors.field_monitors().size() << "\n";
        std::cout << "Harminv probes: " << monitors.harminv_monitors().size() << "\n";
    }

    // Parameter sweep scaffold.
    photonics::SweepRunner sweeper;
    const auto grid = sweeper.make_grid({geom_params.delta_x1}, {geom_params.delta_x2}, {geom_params.edge_radius},
                                        geom_params.lattice_constant);

    if (is_master)
    {
        banner("Sweep Grid");
        std::cout << "Sweep combinations: " << grid.size() << " (dx1/dx2/r_edge)\n";
        for (const auto &p : grid)
        {
            std::cout << "  dx1=" << p.delta_x1 << ", dx2=" << p.delta_x2 << ", r_edge=" << p.r_edge << "\n";
        }
    }

    if (is_master)
    {
        banner("Next Steps");
        std::cout << "- Add flux collection (meep flux boxes) and feed into QAnalyzer for channel breakdown.\n";
        std::cout << "- Drive SweepRunner with a simulation callback that returns SimulationSummary per parameter set.\n";
        std::cout << "- Increase rows/columns and check Q convergence vs resolution/PML/cell size.\n";
    }

    return EXIT_SUCCESS;
}
