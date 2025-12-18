#include "field_slice_export.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace photonics
{
void write_field_slice_csv(const std::string &path, meep::fields &fields, const FieldSliceExportConfig &cfg)
{
    // meep::fields::get_field(..., parallel=true) performs a collective reduction (sum_to_all),
    // so all MPI ranks must call it in the same order. We therefore have every rank execute the
    // same sampling loop, but only rank 0 writes the CSV.
    const bool is_master = !meep::with_mpi() || meep::my_rank() == 0;

    if (cfg.nx <= 1 || cfg.ny <= 1)
    {
        throw std::runtime_error("write_field_slice_csv: nx and ny must be > 1");
    }
    if (!(cfg.sx > 0.0) || !(cfg.sy > 0.0))
    {
        throw std::runtime_error("write_field_slice_csv: sx and sy must be positive");
    }

    std::ofstream out;
    if (is_master)
    {
        out.open(path, std::ios::trunc);
        if (!out)
        {
            throw std::runtime_error("write_field_slice_csv: failed to open for writing: " + path);
        }
    }

    const double xmin = -0.5 * cfg.sx;
    const double xmax = 0.5 * cfg.sx;
    const double ymin = -0.5 * cfg.sy;
    const double ymax = 0.5 * cfg.sy;
    const double dx = cfg.sx / static_cast<double>(cfg.nx);
    const double dy = cfg.sy / static_cast<double>(cfg.ny);
    const double t_now = fields.time();

    if (is_master)
    {
        out << std::setprecision(17);
        out << "# nx=" << cfg.nx
            << " ny=" << cfg.ny
            << " xmin=" << xmin
            << " xmax=" << xmax
            << " ymin=" << ymin
            << " ymax=" << ymax
            << " z=" << cfg.z
            << " component=" << meep::component_name(cfg.component)
            << " part=" << (cfg.imag_part ? "imag" : "real")
            << " time=" << t_now
            << "\n";
        out << std::fixed << std::setprecision(12);
    }
    for (int j = 0; j < cfg.ny; ++j)
    {
        // Sample at pixel centers to stay strictly inside the simulation cell.
        const double y = ymin + (static_cast<double>(j) + 0.5) * dy;

        for (int i = 0; i < cfg.nx; ++i)
        {
            const double x = xmin + (static_cast<double>(i) + 0.5) * dx;
            const auto value = fields.get_field(cfg.component, meep::vec(x, y, cfg.z), true);
            const double scalar = cfg.imag_part ? std::imag(value) : std::real(value);

            if (!is_master)
            {
                continue;
            }

            if (i > 0)
            {
                out << ",";
            }
            out << scalar;
        }
        if (is_master)
        {
            out << "\n";
        }
    }
}

} // namespace photonics
