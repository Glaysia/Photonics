#include "simulation_runner.hpp"

#include <algorithm>
#include <complex>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace photonics
{
namespace
{
double unity_eps(const meep::vec &)
{
    return 1.0;
}

bool is_point_volume(const meep::volume &v)
{
    return v.diameter() <= 1e-12;
}

bool is_master_rank()
{
    return !meep::with_mpi() || meep::my_rank() == 0;
}

std::string format_hms(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0.0)
    {
        seconds = 0.0;
    }
    const long total = static_cast<long>(std::llround(seconds));
    const long h = total / 3600;
    const long m = (total % 3600) / 60;
    const long s = total % 60;

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << ":" << std::setw(2) << m << ":" << std::setw(2) << s;
    return oss.str();
}

} // namespace

SimulationRunResult SimulationRunner::run(const SimulationConfig &config,
                                          const std::vector<meep::geometric_object> &geometry,
                                          const std::vector<meep::source> &sources,
                                          const SimulationRunConfig &run_config) const
{
    SimulationRunResult result{};
    meep::simple_material_function background(&unity_eps);
    meep::grid_volume volume = config.make_grid_volume();
    meep::boundary_region br = config.make_pml();

    result.structure = std::make_unique<meep::structure>(volume, background, br, meep::identity(), 0, config.courant());
    if (!geometry.empty())
    {
        auto owned = geometry;
        meep::geometric_object_list geom_list{};
        geom_list.num_items = static_cast<int>(owned.size());
        geom_list.items = owned.data();
        meep_geom::set_materials_from_geometry(result.structure.get(), geom_list);
    }
    result.fields = std::make_unique<meep::fields>(result.structure.get());

    for (const auto &src : sources)
    {
        if (!src.time)
        {
            throw std::runtime_error("SimulationRunner: source time profile is null");
        }
        if (!is_point_volume(src.where))
        {
            throw std::runtime_error("SimulationRunner: only point sources are supported (volume source not implemented)");
        }

        result.fields->add_point_source(src.component,
                                        *src.time,
                                        src.where.center(),
                                        std::complex<double>(src.amplitude, 0.0));
    }

    result.last_source_time = result.fields->last_source_time();
    result.stop_time = result.last_source_time + run_config.until_after_sources;
    result.steps = 0;
    result.timestep = result.fields->dt;

    const bool collect_ts = run_config.time_series.enabled && run_config.time_series.sample_every_steps >= 1;
    const double ts_start = result.last_source_time + std::max(0.0, run_config.time_series.start_after_sources);
    int ts_stride = run_config.time_series.sample_every_steps;
    if (collect_ts)
    {
        result.sample_dt = result.timestep * static_cast<double>(ts_stride);
    }

    const bool report_progress = is_master_rank() && run_config.progress_interval_seconds > 0.0;
    const double sim_time_start = result.fields->time();
    const double sim_time_total = std::max(result.stop_time - sim_time_start, 0.0);
    const double dt = result.fields->dt;

    using clock = std::chrono::steady_clock;
    const auto wall_start = clock::now();
    auto wall_last = wall_start;
    int steps_last = 0;

    while (result.fields->time() + 0.5 * result.fields->dt < result.stop_time)
    {
        result.fields->step();
        ++result.steps;

        if (collect_ts && result.fields->time() >= ts_start)
        {
            if ((result.steps % ts_stride) == 0)
            {
                const auto value =
                    result.fields->get_field(run_config.time_series.component, run_config.time_series.point, true);
                if (is_master_rank())
                {
                    result.samples.push_back(value);
                }
            }
        }

        if (report_progress)
        {
            const auto wall_now = clock::now();
            const std::chrono::duration<double> wall_since_last = wall_now - wall_last;
            if (wall_since_last.count() >= run_config.progress_interval_seconds)
            {
                const double sim_now = result.fields->time();
                const double sim_done = std::clamp(sim_now - sim_time_start, 0.0, sim_time_total);
                const double progress = sim_time_total > 0.0 ? (sim_done / sim_time_total) : 1.0;

                const int delta_steps = result.steps - steps_last;
                const double steps_per_sec =
                    wall_since_last.count() > 0.0 ? static_cast<double>(delta_steps) / wall_since_last.count() : 0.0;
                const double sec_per_step = steps_per_sec > 0.0 ? 1.0 / steps_per_sec : 0.0;

                const double remaining_time = std::max(result.stop_time - sim_now, 0.0);
                const double remaining_steps = (dt > 0.0) ? std::ceil(remaining_time / dt) : 0.0;
                const double eta_sec = steps_per_sec > 0.0 ? remaining_steps / steps_per_sec : 0.0;

                const double wall_elapsed = std::chrono::duration<double>(wall_now - wall_start).count();

                std::cout << std::fixed << std::setprecision(1);
                std::cout << "[progress] " << (100.0 * progress) << "%  "
                          << "t=" << sim_now << "/" << result.stop_time << "  "
                          << "step=" << result.steps << "  "
                          << steps_per_sec << " steps/s (" << sec_per_step << " s/step)  "
                          << "ETA " << format_hms(eta_sec) << "  "
                          << "elapsed " << format_hms(wall_elapsed) << "\n";
                std::cout.flush();

                wall_last = wall_now;
                steps_last = result.steps;
            }
        }
    }

    if (report_progress)
    {
        const auto wall_end = clock::now();
        const double wall_elapsed = std::chrono::duration<double>(wall_end - wall_start).count();
        std::cout << "[progress] 100.0%  "
                  << "t=" << result.fields->time() << "/" << result.stop_time << "  "
                  << "step=" << result.steps << "  "
                  << "elapsed " << format_hms(wall_elapsed) << "\n";
        std::cout.flush();
    }

    return result;
}

} // namespace photonics
