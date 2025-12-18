
#include "lattice_geometry.hpp"
#include "geometry_export.hpp"
#include "field_slice_export.hpp"
#include "monitor_suite.hpp"
#include "q_analyzer.hpp"
#include "harminv_runner.hpp"
#include "simulation_config.hpp"
#include "simulation_runner.hpp"
#include "source_manager.hpp"
#include "sweep_runner.hpp"

#include <meep.hpp>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
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

std::string json_escape(const std::string &s)
{
    std::ostringstream oss;
    for (const char c : s)
    {
        switch (c)
        {
        case '\\':
            oss << "\\\\";
            break;
        case '"':
            oss << "\\\"";
            break;
        case '\n':
            oss << "\\n";
            break;
        case '\r':
            oss << "\\r";
            break;
        case '\t':
            oss << "\\t";
            break;
        default:
            oss << c;
        }
    }
    return oss.str();
}

std::string iso8601_utc_now()
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &tt);
#else
    tm = *std::gmtime(&tt);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

std::string lower_ascii(std::string s)
{
    for (char &c : s)
    {
        if (c >= 'A' && c <= 'Z')
        {
            c = static_cast<char>(c - 'A' + 'a');
        }
    }
    return s;
}

bool parse_bool(const std::string &v, bool fallback)
{
    const std::string s = lower_ascii(v);
    if (s == "1" || s == "true" || s == "yes" || s == "on")
    {
        return true;
    }
    if (s == "0" || s == "false" || s == "no" || s == "off")
    {
        return false;
    }
    return fallback;
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

int process_size()
{
    if (meep::with_mpi())
    {
        return meep::count_processors();
    }
    for (const char *name : {"OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "SLURM_NTASKS", "MV2_COMM_WORLD_SIZE"})
    {
        const int v = parse_int_env(name);
        if (v > 0)
        {
            return v;
        }
    }
    return 1;
}

int omp_threads()
{
    const int v = parse_int_env("OMP_NUM_THREADS");
    if (v > 0)
    {
        return v;
    }
    return 1;
}

std::string format_time_tag(double t)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << t;
    std::string s = oss.str();
    while (!s.empty() && s.back() == '0')
    {
        s.pop_back();
    }
    if (!s.empty() && s.back() == '.')
    {
        s.pop_back();
    }
    if (!s.empty() && s.front() == '-')
    {
        s = std::string("m") + s.substr(1);
    }
    for (char &c : s)
    {
        if (c == '.')
        {
            c = 'p';
        }
    }
    return s.empty() ? std::string("0") : s;
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

struct DumpProbe
{
    std::string label;
    meep::vec point{0.0, 0.0, 0.0};
    meep::component component = meep::Ey;
};

struct DumpTimeSeries
{
    int sample_every_steps = 0;
    double dt = 0.0;
    double sample_dt = 0.0;
    double start_time = 0.0;
    std::vector<double> times;
    std::vector<std::vector<std::complex<double>>> values; // [probe][sample]
};

struct DumpFluxSpectrum
{
    double fmin = 0.0;
    double fmax = 0.0;
    int nf = 0;
    std::vector<double> freq;
    std::vector<double> top;
    std::vector<double> bottom;
    std::vector<double> left;
    std::vector<double> right;
    std::vector<double> front;
    std::vector<double> back;
    std::vector<double> total;
};

struct DumpSimulationResult
{
    double timestep = 0.0;
    double last_source_time = 0.0;
    double stop_time = 0.0;
    int steps = 0;

    DumpTimeSeries time_series;
    DumpFluxSpectrum flux;

    // Saved (requested) snapshot outputs: (target_time, path_written)
    std::vector<std::pair<double, std::string>> snapshots;
};

bool is_point_volume(const meep::volume &v)
{
    return v.diameter() <= 1e-12;
}

std::string component_label(meep::component c)
{
    return meep::component_name(c);
}

void write_time_series_csv(const std::string &path, const DumpSimulationResult &res, const std::vector<DumpProbe> &probes)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("failed to open for writing: " + path);
    }

    out << std::setprecision(17);
    out << "# dt=" << res.time_series.dt
        << " sample_every_steps=" << res.time_series.sample_every_steps
        << " sample_dt=" << res.time_series.sample_dt
        << " start_time=" << res.time_series.start_time
        << "\n";
    for (std::size_t i = 0; i < probes.size(); ++i)
    {
        const auto &p = probes[i];
        out << "# probe" << i << " label=" << p.label
            << " component=" << component_label(p.component)
            << " x=" << p.point.x()
            << " y=" << p.point.y()
            << " z=" << p.point.z()
            << "\n";
    }

    out << "t";
    for (std::size_t i = 0; i < probes.size(); ++i)
    {
        out << ",re" << i << ",im" << i;
    }
    out << "\n";

    const std::size_t n = res.time_series.times.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        out << res.time_series.times[k];
        for (std::size_t i = 0; i < probes.size(); ++i)
        {
            const auto &v = res.time_series.values[i][k];
            out << "," << std::real(v) << "," << std::imag(v);
        }
        out << "\n";
    }
}

void write_flux_spectrum_csv(const std::string &path, const DumpFluxSpectrum &flux)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("failed to open for writing: " + path);
    }

    out << std::setprecision(17);
    out << "# fmin=" << flux.fmin << " fmax=" << flux.fmax << " nf=" << flux.nf << "\n";
    out << "freq,top,bottom,left,right,front,back,total\n";
    const std::size_t n = flux.freq.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << flux.freq[i] << ","
            << flux.top[i] << ","
            << flux.bottom[i] << ","
            << flux.left[i] << ","
            << flux.right[i] << ","
            << flux.front[i] << ","
            << flux.back[i] << ","
            << flux.total[i] << "\n";
    }
}

DumpSimulationResult run_dump_time_domain(const photonics::SimulationConfig &config,
                                         const std::vector<meep::geometric_object> &geometry,
                                         const std::vector<meep::source> &sources,
                                         const photonics::SimulationParameters &sim_params,
                                         const photonics::SourceConfig &source_cfg,
                                         const std::tuple<int, int, int> &grid_points,
                                         const std::string &dump_prefix,
                                         double until_after_sources,
                                         double progress_interval_seconds,
                                         const std::vector<DumpProbe> &probes,
                                         int ts_stride,
                                         double ts_start_after_sources)
{
    DumpSimulationResult result{};

    // Build structure/fields.
    meep::structure structure = config.make_structure(geometry);
    meep::fields fields(&structure);

    // Add point sources.
    for (const auto &src : sources)
    {
        if (!src.time)
        {
            throw std::runtime_error("dump: source time profile is null");
        }
        if (!is_point_volume(src.where))
        {
            throw std::runtime_error("dump: only point sources are supported");
        }
        fields.add_point_source(src.component,
                                *src.time,
                                src.where.center(),
                                std::complex<double>(src.amplitude, 0.0));
    }

    result.last_source_time = fields.last_source_time();
    result.stop_time = result.last_source_time + until_after_sources;
    result.timestep = fields.dt;
    result.steps = 0;

    // Configure time-series storage (master rank only stores the samples).
    result.time_series.sample_every_steps = std::max(1, ts_stride);
    result.time_series.dt = result.timestep;
    result.time_series.sample_dt = result.timestep * static_cast<double>(result.time_series.sample_every_steps);
    result.time_series.start_time = result.last_source_time + std::max(0.0, ts_start_after_sources);
    result.time_series.values.resize(probes.size());

    // Snapshot exports (all ranks must call in same order; only rank 0 writes).
    photonics::FieldSliceExportConfig slice_cfg{};
    slice_cfg.component = meep::Ey;
    slice_cfg.z = 0.0;
    slice_cfg.nx = std::max(32, std::get<0>(grid_points));
    slice_cfg.ny = std::max(32, std::get<1>(grid_points));
    slice_cfg.sx = config.cell_size_x_in_a();
    slice_cfg.sy = config.cell_size_y_in_a();
    slice_cfg.imag_part = false;

    std::vector<double> snapshot_times;
    snapshot_times.reserve(3);
    if (source_cfg.bandwidth > 0.0)
    {
        const double peak_time = source_cfg.cutoff / source_cfg.bandwidth;
        snapshot_times.push_back(peak_time);
    }
    snapshot_times.push_back(result.last_source_time);
    snapshot_times.push_back(result.last_source_time + 0.5 * until_after_sources);
    // Deduplicate (within half a timestep).
    std::vector<double> snapshots;
    for (const double t : snapshot_times)
    {
        if (!(t >= 0.0) || !(t <= result.stop_time))
        {
            continue;
        }
        bool dup = false;
        for (const double existing : snapshots)
        {
            if (std::abs(existing - t) <= 0.5 * result.timestep)
            {
                dup = true;
                break;
            }
        }
        if (!dup)
        {
            snapshots.push_back(t);
        }
    }
    std::vector<bool> snapshot_done(snapshots.size(), false);

    // DFT flux monitors.
    const double fc = source_cfg.center_frequency;
    const double bw = std::max(source_cfg.bandwidth, 1e-6);
    result.flux.fmin = fc - 2.0 * bw;
    result.flux.fmax = fc + 2.0 * bw;
    result.flux.nf = 200;
    result.flux.freq = meep::linspace(result.flux.fmin, result.flux.fmax, static_cast<std::size_t>(result.flux.nf));

    photonics::MonitorConfig monitor_cfg{};
    monitor_cfg.cell = {sim_params.cell_size_x_in_a, sim_params.cell_size_y_in_a, sim_params.cell_size_z_in_a, meep::D3};
    monitor_cfg.pml_thickness = sim_params.pml_thickness_in_a;
    monitor_cfg.monitor_thickness = 0.05;
    photonics::MonitorSuite monitors(monitor_cfg);
    monitors.setup_default_flux_boxes();

    auto flux_region = [&](const std::string &label) -> meep::volume {
        for (const auto &m : monitors.flux_monitors())
        {
            if (m.label == label)
            {
                return m.region;
            }
        }
        throw std::runtime_error("dump: missing flux monitor region: " + label);
    };

    meep::dft_flux flux_left = fields.add_dft_flux(meep::X, flux_region("flux-left"), result.flux.freq);
    meep::dft_flux flux_right = fields.add_dft_flux(meep::X, flux_region("flux-right"), result.flux.freq);
    meep::dft_flux flux_front = fields.add_dft_flux(meep::Y, flux_region("flux-front"), result.flux.freq);
    meep::dft_flux flux_back = fields.add_dft_flux(meep::Y, flux_region("flux-back"), result.flux.freq);
    meep::dft_flux flux_top = fields.add_dft_flux(meep::Z, flux_region("flux-top"), result.flux.freq);
    meep::dft_flux flux_bottom = fields.add_dft_flux(meep::Z, flux_region("flux-bottom"), result.flux.freq);

    const bool is_master = (process_rank() == 0);
    using clock = std::chrono::steady_clock;
    const auto wall_start = clock::now();
    auto wall_last = wall_start;
    int steps_last = 0;

    const double sim_time_start = fields.time();
    const double sim_time_total = std::max(result.stop_time - sim_time_start, 0.0);

    while (fields.time() + 0.5 * fields.dt < result.stop_time)
    {
        fields.step();
        ++result.steps;

        // Time series sampling (collective get_field; only rank 0 stores values).
        if (fields.time() >= result.time_series.start_time && (result.steps % result.time_series.sample_every_steps) == 0)
        {
            if (is_master)
            {
                result.time_series.times.push_back(fields.time());
            }
            for (std::size_t i = 0; i < probes.size(); ++i)
            {
                const auto v = fields.get_field(probes[i].component, probes[i].point, true);
                if (is_master)
                {
                    result.time_series.values[i].push_back(v);
                }
            }
        }

        // Snapshots (collective get_field inside write_field_slice_csv).
        for (std::size_t i = 0; i < snapshots.size(); ++i)
        {
            if (snapshot_done[i])
            {
                continue;
            }
            if (fields.time() + 0.5 * fields.dt < snapshots[i])
            {
                continue;
            }
            const std::string out_path = dump_prefix + "_Ey_z0_t" + format_time_tag(snapshots[i]) + ".csv";
            if (is_master)
            {
                ensure_parent_dir(out_path);
            }
            photonics::write_field_slice_csv(out_path, fields, slice_cfg);
            result.snapshots.emplace_back(snapshots[i], out_path);
            snapshot_done[i] = true;
        }

        // Periodic progress logging.
        if (is_master && progress_interval_seconds > 0.0)
        {
            const auto wall_now = clock::now();
            const std::chrono::duration<double> wall_since_last = wall_now - wall_last;
            if (wall_since_last.count() >= progress_interval_seconds)
            {
                const double sim_now = fields.time();
                const double sim_done = std::clamp(sim_now - sim_time_start, 0.0, sim_time_total);
                const double progress = sim_time_total > 0.0 ? (sim_done / sim_time_total) : 1.0;

                const int delta_steps = result.steps - steps_last;
                const double steps_per_sec =
                    wall_since_last.count() > 0.0 ? static_cast<double>(delta_steps) / wall_since_last.count() : 0.0;
                const double sec_per_step = steps_per_sec > 0.0 ? 1.0 / steps_per_sec : 0.0;

                const double remaining_time = std::max(result.stop_time - sim_now, 0.0);
                const double remaining_steps = (fields.dt > 0.0) ? std::ceil(remaining_time / fields.dt) : 0.0;
                const double eta_sec = steps_per_sec > 0.0 ? remaining_steps / steps_per_sec : 0.0;

                const double wall_elapsed = std::chrono::duration<double>(wall_now - wall_start).count();

                std::cout << std::fixed << std::setprecision(1);
                std::cout << "[dump progress] " << (100.0 * progress) << "%  "
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

    // Collect DFT fluxes (collective).
    auto collect_flux = [&](meep::dft_flux &f) -> std::vector<double> {
        double *arr = f.flux();
        std::vector<double> v;
        v.reserve(f.freq.size());
        for (std::size_t i = 0; i < f.freq.size(); ++i)
        {
            v.push_back(arr[i]);
        }
        delete[] arr;
        return v;
    };

    // Raw face fluxes in the +normal direction (Meep convention).
    const auto f_left = collect_flux(flux_left);
    const auto f_right = collect_flux(flux_right);
    const auto f_front = collect_flux(flux_front);
    const auto f_back = collect_flux(flux_back);
    const auto f_top = collect_flux(flux_top);
    const auto f_bottom = collect_flux(flux_bottom);

    // Convert to outward flux for a closed box.
    const std::size_t nfreq = result.flux.freq.size();
    result.flux.left.resize(nfreq);
    result.flux.right.resize(nfreq);
    result.flux.front.resize(nfreq);
    result.flux.back.resize(nfreq);
    result.flux.top.resize(nfreq);
    result.flux.bottom.resize(nfreq);
    result.flux.total.resize(nfreq);
    for (std::size_t i = 0; i < nfreq; ++i)
    {
        result.flux.right[i] = f_right[i];
        result.flux.left[i] = -f_left[i];
        result.flux.front[i] = f_front[i];
        result.flux.back[i] = -f_back[i];
        result.flux.top[i] = f_top[i];
        result.flux.bottom[i] = -f_bottom[i];
        result.flux.total[i] = result.flux.top[i] + result.flux.bottom[i] + result.flux.left[i] + result.flux.right[i] +
                               result.flux.front[i] + result.flux.back[i];
    }

    return result;
}

struct SelectedHarminvMode
{
    bool used_fallback = false;
    std::string method;
    std::size_t index = 0;
    double frequency = 0.0;
};

SelectedHarminvMode select_harminv_mode(const std::vector<photonics::HarminvMode> &modes, double fallback_frequency)
{
    SelectedHarminvMode selected{};
    if (modes.empty())
    {
        selected.used_fallback = true;
        selected.method = "fallback_center_frequency";
        selected.frequency = fallback_frequency;
        return selected;
    }

    double best_q = -1.0;
    double best_amp = -1.0;
    std::size_t best = 0;
    bool found = false;
    for (std::size_t i = 0; i < modes.size(); ++i)
    {
        const double freq = modes[i].frequency;
        const double decay = std::abs(modes[i].decay_rate);
        if (!(freq > 0.0) || !(decay > 0.0))
        {
            continue;
        }
        const double q = freq / (2.0 * decay);
        const double amp = std::abs(modes[i].amplitude);
        if (!found || q > best_q || (std::abs(q - best_q) < 1e-12 && amp > best_amp))
        {
            found = true;
            best_q = q;
            best_amp = amp;
            best = i;
        }
    }
    if (!found)
    {
        selected.used_fallback = true;
        selected.method = "fallback_center_frequency";
        selected.frequency = fallback_frequency;
        return selected;
    }

    selected.used_fallback = false;
    selected.method = "max_Q";
    selected.index = best;
    selected.frequency = modes[best].frequency;
    return selected;
}

void write_harminv_json(const std::string &path,
                        const DumpProbe &probe,
                        const DumpSimulationResult &run,
                        const photonics::HarminvRunConfig &config,
                        const std::vector<photonics::HarminvMode> &modes,
                        const SelectedHarminvMode &selected)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("failed to open for writing: " + path);
    }

    out << std::fixed << std::setprecision(12);
    out << "{\n";
    out << "  \"probe\": {\n";
    out << "    \"label\": \"" << json_escape(probe.label) << "\",\n";
    out << "    \"component\": \"" << json_escape(component_label(probe.component)) << "\",\n";
    out << "    \"point\": [" << probe.point.x() << ", " << probe.point.y() << ", " << probe.point.z() << "],\n";
    out << "    \"dt\": " << run.time_series.dt << ",\n";
    out << "    \"sample_every_steps\": " << run.time_series.sample_every_steps << ",\n";
    out << "    \"sample_dt\": " << run.time_series.sample_dt << ",\n";
    out << "    \"start_time\": " << run.time_series.start_time << ",\n";
    out << "    \"num_samples\": " << run.time_series.times.size() << "\n";
    out << "  },\n";
    out << "  \"config\": {\n";
    out << "    \"fmin\": " << config.fmin << ",\n";
    out << "    \"fmax\": " << config.fmax << ",\n";
    out << "    \"maxbands\": " << config.maxbands << ",\n";
    out << "    \"q_threshold\": " << config.q_threshold << "\n";
    out << "  },\n";
    out << "  \"selected\": {\n";
    out << "    \"method\": \"" << json_escape(selected.method) << "\",\n";
    if (selected.used_fallback)
    {
        out << "    \"used_fallback\": true,\n";
        out << "    \"frequency\": " << selected.frequency << "\n";
    }
    else
    {
        const auto &m = modes[selected.index];
        const double decay = std::abs(m.decay_rate);
        const double q = decay > 0.0 ? m.frequency / (2.0 * decay) : 0.0;
        out << "    \"used_fallback\": false,\n";
        out << "    \"index\": " << selected.index << ",\n";
        out << "    \"frequency\": " << m.frequency << ",\n";
        out << "    \"decay\": " << decay << ",\n";
        out << "    \"Q\": " << q << ",\n";
        out << "    \"amp\": " << std::abs(m.amplitude) << ",\n";
        out << "    \"phase\": " << m.phase << "\n";
    }
    out << "  },\n";
    out << "  \"modes\": [\n";
    for (std::size_t i = 0; i < modes.size(); ++i)
    {
        const auto &m = modes[i];
        const double decay = std::abs(m.decay_rate);
        const double q = decay > 0.0 ? m.frequency / (2.0 * decay) : 0.0;
        out << "    {\"freq\": " << m.frequency
            << ", \"decay\": " << decay
            << ", \"Q\": " << q
            << ", \"amp\": " << std::abs(m.amplitude)
            << ", \"phase\": " << m.phase << "}";
        if (i + 1 != modes.size())
        {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void dump_dft_ey_z0_csv(const std::string &path,
                        const photonics::SimulationConfig &config,
                        const std::vector<meep::geometric_object> &geometry,
                        const std::vector<meep::source> &sources,
                        double frequency,
                        double after_sources_time)
{
    meep::structure structure = config.make_structure(geometry);
    meep::fields fields(&structure);

    for (const auto &src : sources)
    {
        if (!src.time)
        {
            throw std::runtime_error("dft: source time profile is null");
        }
        if (!is_point_volume(src.where))
        {
            throw std::runtime_error("dft: only point sources are supported");
        }
        fields.add_point_source(src.component,
                                *src.time,
                                src.where.center(),
                                std::complex<double>(src.amplitude, 0.0));
    }

    const double last_source_time = fields.last_source_time();
    const double stop_time = last_source_time + std::max(0.0, after_sources_time);

    const double sx = config.cell_size_x_in_a();
    const double sy = config.cell_size_y_in_a();
    const double z = 0.0;
    const meep::volume plane(meep::vec(-0.5 * sx, -0.5 * sy, z), meep::vec(0.5 * sx, 0.5 * sy, z));
    const std::vector<double> freqs{frequency};
    meep::component comps[1] = {meep::Ey};

    // Excite the cavity with the broadband source first.
    while (fields.time() + 0.5 * fields.dt < last_source_time)
    {
        fields.step();
    }

    // Start the frequency-domain accumulation after sources have ended (ringdown only).
    meep::dft_fields dft = fields.add_dft_fields(comps, 1, plane, freqs);

    while (fields.time() + 0.5 * fields.dt < stop_time)
    {
        fields.step();
    }

    int rank = 0;
    size_t dims[3] = {0, 0, 0};
    auto *arr = fields.get_dft_array(dft, meep::Ey, 0, &rank, dims);

    const bool is_master = (process_rank() == 0);
    if (is_master)
    {
        ensure_parent_dir(path);
        std::ofstream out(path, std::ios::trunc);
        if (!out)
        {
            delete[] arr;
            throw std::runtime_error("failed to open for writing: " + path);
        }

        out << std::setprecision(17);
        out << "# nx=" << dims[0]
            << " ny=" << (rank >= 2 ? dims[1] : 1)
            << " xmin=" << (-0.5 * sx)
            << " xmax=" << (0.5 * sx)
            << " ymin=" << (-0.5 * sy)
            << " ymax=" << (0.5 * sy)
            << " z=" << z
            << " component=" << component_label(meep::Ey)
            << " freq=" << frequency
            << "\n";
        out << "x,y,re,im,abs\n";

        const std::size_t nx = dims[0];
        const std::size_t ny = (rank >= 2 ? dims[1] : 1);
        const double xmin = -0.5 * sx;
        const double ymin = -0.5 * sy;
        const double dx = sx / static_cast<double>(nx);
        const double dy = sy / static_cast<double>(ny);
        for (std::size_t j = 0; j < ny; ++j)
        {
            const double yv = ymin + (static_cast<double>(j) + 0.5) * dy;
            for (std::size_t i = 0; i < nx; ++i)
            {
                const double xv = xmin + (static_cast<double>(i) + 0.5) * dx;
                const auto v = arr[i * ny + j];
                out << xv << "," << yv << "," << std::real(v) << "," << std::imag(v) << "," << std::abs(v) << "\n";
            }
        }
    }

    delete[] arr;
}

void write_run_json(const std::string &path,
                    int argc,
                    char **argv,
                    const std::string &variant,
                    const std::string &dump_prefix,
                    const photonics::SimulationParameters &sim_params,
                    const photonics::SourceConfig &source_cfg,
                    const DumpSimulationResult &run,
                    const photonics::HarminvRunConfig &hcfg,
                    const SelectedHarminvMode &selected,
                    double dft_after_sources_time)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("failed to open for writing: " + path);
    }

    out << std::fixed << std::setprecision(12);
    out << "{\n";
    out << "  \"timestamp_utc\": \"" << iso8601_utc_now() << "\",\n";
    out << "  \"argv\": [";
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << "\"" << json_escape(argv[i]) << "\"";
    }
    out << "],\n";
    out << "  \"mpi\": {\"np\": " << process_size() << ", \"rank\": " << process_rank() << "},\n";
    out << "  \"omp\": {\"num_threads\": " << omp_threads() << "},\n";
    out << "  \"variant\": \"" << json_escape(variant) << "\",\n";
    out << "  \"dump_prefix\": \"" << json_escape(dump_prefix) << "\",\n";

    out << "  \"simulation_parameters\": {\n";
    out << "    \"lattice_constant_um\": " << sim_params.lattice_constant_um << ",\n";
    out << "    \"resolution_px_per_a\": " << sim_params.resolution_px_per_a << ",\n";
    out << "    \"cell_size_x_in_a\": " << sim_params.cell_size_x_in_a << ",\n";
    out << "    \"cell_size_y_in_a\": " << sim_params.cell_size_y_in_a << ",\n";
    out << "    \"cell_size_z_in_a\": " << sim_params.cell_size_z_in_a << ",\n";
    out << "    \"pml_thickness_in_a\": " << sim_params.pml_thickness_in_a << ",\n";
    out << "    \"courant\": " << sim_params.courant << "\n";
    out << "  },\n";

    out << "  \"source_config\": {\n";
    out << "    \"center_frequency\": " << source_cfg.center_frequency << ",\n";
    out << "    \"bandwidth\": " << source_cfg.bandwidth << ",\n";
    out << "    \"cutoff\": " << source_cfg.cutoff << ",\n";
    out << "    \"amplitude\": " << source_cfg.amplitude << ",\n";
    out << "    \"components\": [";
    for (std::size_t i = 0; i < source_cfg.components.size(); ++i)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << "\"" << json_escape(component_label(source_cfg.components[i])) << "\"";
    }
    out << "],\n";
    out << "    \"positions\": [";
    for (std::size_t i = 0; i < source_cfg.positions.size(); ++i)
    {
        const auto &p = source_cfg.positions[i];
        if (i > 0)
        {
            out << ", ";
        }
        out << "[" << p.x() << ", " << p.y() << ", " << p.z() << "]";
    }
    out << "]\n";
    out << "  },\n";

    out << "  \"run\": {\n";
    out << "    \"dt\": " << run.timestep << ",\n";
    out << "    \"steps\": " << run.steps << ",\n";
    out << "    \"last_source_time\": " << run.last_source_time << ",\n";
    out << "    \"stop_time\": " << run.stop_time << "\n";
    out << "  },\n";

    out << "  \"time_series\": {\n";
    out << "    \"sample_every_steps\": " << run.time_series.sample_every_steps << ",\n";
    out << "    \"sample_dt\": " << run.time_series.sample_dt << ",\n";
    out << "    \"start_time\": " << run.time_series.start_time << ",\n";
    out << "    \"num_samples\": " << run.time_series.times.size() << "\n";
    out << "  },\n";

    out << "  \"harminv\": {\n";
    out << "    \"fmin\": " << hcfg.fmin << ",\n";
    out << "    \"fmax\": " << hcfg.fmax << ",\n";
    out << "    \"maxbands\": " << hcfg.maxbands << ",\n";
    out << "    \"q_threshold\": " << hcfg.q_threshold << ",\n";
    out << "    \"selected_frequency\": " << selected.frequency << ",\n";
    out << "    \"selection_method\": \"" << json_escape(selected.method) << "\"\n";
    out << "  },\n";

    out << "  \"dft_field_map\": {\n";
    out << "    \"component\": \"" << json_escape(component_label(meep::Ey)) << "\",\n";
    out << "    \"plane\": {\"normal\": \"Z\", \"coordinate\": 0.0},\n";
    out << "    \"frequency\": " << selected.frequency << ",\n";
    out << "    \"after_sources_time\": " << dft_after_sources_time << "\n";
    out << "  },\n";

    out << "  \"flux_spectrum\": {\n";
    out << "    \"fmin\": " << run.flux.fmin << ",\n";
    out << "    \"fmax\": " << run.flux.fmax << ",\n";
    out << "    \"nf\": " << run.flux.nf << "\n";
    out << "  },\n";

    out << "  \"snapshots\": [\n";
    for (std::size_t i = 0; i < run.snapshots.size(); ++i)
    {
        out << "    {\"target_time\": " << run.snapshots[i].first << ", \"path\": \"" << json_escape(run.snapshots[i].second)
            << "\"}";
        if (i + 1 != run.snapshots.size())
        {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
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

    const std::string variant = [&]() {
        if (has_flag(argc, argv, "--baseline") || has_flag(argc, argv, "--unshifted"))
        {
            return std::string("baseline");
        }
        if (has_flag(argc, argv, "--shifted"))
        {
            return std::string("shifted");
        }
        const auto v = lower_ascii(arg_value(argc, argv, "--variant"));
        if (v == "baseline" || v == "unshifted" || v == "0")
        {
            return std::string("baseline");
        }
        if (v == "shifted" || v == "1")
        {
            return std::string("shifted");
        }
        return std::string("shifted");
    }();

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
    geom_params.delta_x1 = (variant == "shifted") ? (0.15 * geom_params.lattice_constant) : 0.0;
    geom_params.delta_x2 = (variant == "shifted") ? (0.05 * geom_params.lattice_constant) : 0.0;
    geom_params.edge_radius = geom_params.hole_radius;
    if (const auto v = arg_value(argc, argv, "--edge-radius"); !v.empty())
    {
        geom_params.edge_radius = std::stod(v);
    }

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
        meta.units = "a";

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

    const std::string dump_prefix = [&]() {
        // New preferred flag: `--dump <prefix>`.
        if (const auto v = arg_value(argc, argv, "--dump"); !v.empty())
        {
            return v;
        }
        // Backward-compatible alias: `--dump-colormap [prefix]` maps to `--dump` behavior.
        const auto v = arg_value(argc, argv, "--dump-colormap");
        if (!v.empty())
        {
            return v;
        }
        if (has_flag(argc, argv, "--dump-colormap"))
        {
            return std::string("out/") + variant;
        }
        return std::string();
    }();

    if (has_flag(argc, argv, "--dump") && dump_prefix.empty())
    {
        if (is_master)
        {
            std::cerr << "error: --dump requires an output prefix, e.g. --dump out/run1\n";
        }
        return EXIT_FAILURE;
    }

    if (is_master && !dump_prefix.empty())
    {
        photonics::GeometryMetadata meta{};
        meta.params = geom_params;
        meta.holes = holes;
        meta.units = "a";

        const auto csv_path = ensure_parent_dir(dump_prefix + ".csv");
        const auto json_path = ensure_parent_dir(dump_prefix + ".json");
        const auto svg_path = ensure_parent_dir(dump_prefix + ".svg");

        photonics::write_holes_csv(csv_path.string(), holes);
        photonics::write_geometry_json(json_path.string(), meta);
        photonics::write_holes_svg(svg_path.string(), holes);
        std::cout << "Exported geometry:\n";
        std::cout << "  - " << csv_path << "\n";
        std::cout << "  - " << json_path << "\n";
        std::cout << "  - " << svg_path << "\n";
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

    if (!dump_prefix.empty())
    {
        try
        {
            const double default_until_after_sources = 200.0;
            double until_after_sources = default_until_after_sources;
            if (const auto v = arg_value(argc, argv, "--run-until-after-sources"); !v.empty())
            {
                until_after_sources = std::stod(v);
            }

            const std::vector<DumpProbe> probes{
                {"center", meep::vec(0.0, 0.0, 0.0), meep::Ey},
                {"x+0.25", meep::vec(0.25, 0.0, 0.0), meep::Ey},
                {"y+0.25", meep::vec(0.0, 0.25, 0.0), meep::Ey},
            };

            if (is_master)
            {
                banner("Dump");
            }

            const DumpSimulationResult dump_run =
                run_dump_time_domain(sim_config,
                                     geometry,
                                     sources,
                                     sim_params,
                                     source_manager.config(),
                                     grid_points,
                                     dump_prefix,
                                     until_after_sources,
                                     /*progress_interval_seconds=*/7.0,
                                     probes,
                                     /*ts_stride=*/10,
                                     /*ts_start_after_sources=*/0.0);

            if (is_master)
            {
                banner("Dump Outputs");
                ensure_parent_dir(dump_prefix + "_ts.csv");
                ensure_parent_dir(dump_prefix + "_flux_spectrum.csv");
                ensure_parent_dir(dump_prefix + "_harminv.json");
                ensure_parent_dir(dump_prefix + "_run.json");

                write_time_series_csv(dump_prefix + "_ts.csv", dump_run, probes);
                write_flux_spectrum_csv(dump_prefix + "_flux_spectrum.csv", dump_run.flux);
            }

            photonics::HarminvRunConfig hcfg{};
            const double fc = source_manager.config().center_frequency;
            const double bw = std::max(source_manager.config().bandwidth, 1e-6);
            hcfg.fmin = fc - 2.0 * bw;
            hcfg.fmax = fc + 2.0 * bw;
            hcfg.maxbands = 50;
            hcfg.q_threshold = 100.0;

            std::vector<photonics::HarminvMode> modes;
            SelectedHarminvMode selected{};
            if (is_master)
            {
                modes = photonics::run_harminv(dump_run.time_series.values.front(), dump_run.time_series.sample_dt, hcfg);
                selected = select_harminv_mode(modes, fc);
            }

            double selected_freq = is_master ? selected.frequency : 0.0;
            if (meep::with_mpi())
            {
                selected_freq = meep::broadcast(0, selected_freq);
            }
            if (is_master)
            {
                selected.frequency = selected_freq;
            }

            if (is_master)
            {
                write_harminv_json(dump_prefix + "_harminv.json", probes.front(), dump_run, hcfg, modes, selected);
            }

            const double dft_after_sources_time = std::min(until_after_sources, 50.0);
            dump_dft_ey_z0_csv(dump_prefix + "_Ey_z0_dft.csv",
                               sim_config,
                               geometry,
                               sources,
                               selected_freq,
                               dft_after_sources_time);

            if (is_master)
            {
                write_run_json(dump_prefix + "_run.json",
                               argc,
                               argv,
                               variant,
                               dump_prefix,
                               sim_params,
                               source_manager.config(),
                               dump_run,
                               hcfg,
                               selected,
                               dft_after_sources_time);

                std::cout << "Wrote dump prefix: " << dump_prefix << "\n";
            }

            return EXIT_SUCCESS;
        }
        catch (const std::exception &e)
        {
            meep::abort("dump failed: %s", e.what());
        }
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
