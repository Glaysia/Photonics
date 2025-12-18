#include "field_slice_export.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace photonics
{
namespace
{
double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}
} // namespace

void write_field_slice_csv(const std::string &path, meep::fields &fields, const FieldSliceExportConfig &cfg)
{
    if (cfg.nx <= 1 || cfg.ny <= 1)
    {
        throw std::runtime_error("write_field_slice_csv: nx and ny must be > 1");
    }
    if (!(cfg.sx > 0.0) || !(cfg.sy > 0.0))
    {
        throw std::runtime_error("write_field_slice_csv: sx and sy must be positive");
    }

    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        throw std::runtime_error("write_field_slice_csv: failed to open for writing: " + path);
    }

    const double xmin = -0.5 * cfg.sx;
    const double xmax = 0.5 * cfg.sx;
    const double ymin = -0.5 * cfg.sy;
    const double ymax = 0.5 * cfg.sy;
    const double t_now = fields.time();

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
    for (int j = 0; j < cfg.ny; ++j)
    {
        const double ty = static_cast<double>(j) / static_cast<double>(cfg.ny - 1);
        const double y = lerp(ymin, ymax, ty);

        for (int i = 0; i < cfg.nx; ++i)
        {
            const double tx = static_cast<double>(i) / static_cast<double>(cfg.nx - 1);
            const double x = lerp(xmin, xmax, tx);
            const auto value = fields.get_field(cfg.component, meep::vec(x, y, cfg.z), true);
            const double scalar = cfg.imag_part ? std::imag(value) : std::real(value);

            if (i > 0)
            {
                out << ",";
            }
            out << scalar;
        }
        out << "\n";
    }
}

} // namespace photonics
