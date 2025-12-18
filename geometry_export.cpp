#include "geometry_export.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace photonics
{
namespace
{
struct Bounds
{
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
};

Bounds compute_bounds(const std::vector<Hole> &holes)
{
    Bounds b{};
    for (const auto &h : holes)
    {
        b.min_x = std::min(b.min_x, h.center.x() - h.radius);
        b.max_x = std::max(b.max_x, h.center.x() + h.radius);
        b.min_y = std::min(b.min_y, h.center.y() - h.radius);
        b.max_y = std::max(b.max_y, h.center.y() + h.radius);
    }
    if (!std::isfinite(b.min_x) || !std::isfinite(b.min_y))
    {
        b.min_x = b.min_y = -1.0;
        b.max_x = b.max_y = 1.0;
    }
    return b;
}

double safe_div(double num, double den)
{
    return (std::abs(den) < 1e-12) ? 0.0 : (num / den);
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

void require_stream(std::ofstream &out, const std::string &path)
{
    if (!out)
    {
        throw std::runtime_error("Failed to open for writing: " + path);
    }
}

} // namespace

void write_holes_csv(const std::string &path, const std::vector<Hole> &holes)
{
    std::ofstream out(path, std::ios::trunc);
    require_stream(out, path);

    out << "x,y,r\n";
    out << std::fixed << std::setprecision(9);
    for (const auto &h : holes)
    {
        out << h.center.x() << "," << h.center.y() << "," << h.radius << "\n";
    }
}

void write_geometry_json(const std::string &path, const GeometryMetadata &metadata)
{
    std::ofstream out(path, std::ios::trunc);
    require_stream(out, path);

    out << std::fixed << std::setprecision(12);
    out << "{\n";
    out << "  \"units\": \"" << json_escape(metadata.units) << "\",\n";
    out << "  \"params\": {\n";
    out << "    \"lattice_constant\": " << metadata.params.lattice_constant << ",\n";
    out << "    \"hole_radius\": " << metadata.params.hole_radius << ",\n";
    out << "    \"slab_thickness\": " << metadata.params.slab_thickness << ",\n";
    out << "    \"delta_x1\": " << metadata.params.delta_x1 << ",\n";
    out << "    \"delta_x2\": " << metadata.params.delta_x2 << ",\n";
    out << "    \"edge_radius\": " << metadata.params.edge_radius << ",\n";
    out << "    \"slab_index\": " << metadata.params.slab_index << ",\n";
    out << "    \"hole_index\": " << metadata.params.hole_index << ",\n";
    out << "    \"columns\": " << metadata.params.columns << ",\n";
    out << "    \"rows\": " << metadata.params.rows << "\n";
    out << "  },\n";
    out << "  \"holes\": [\n";

    for (std::size_t i = 0; i < metadata.holes.size(); ++i)
    {
        const auto &h = metadata.holes[i];
        out << "    {\"x\": " << h.center.x() << ", \"y\": " << h.center.y() << ", \"r\": " << h.radius << "}";
        if (i + 1 != metadata.holes.size())
        {
            out << ",";
        }
        out << "\n";
    }

    out << "  ]\n";
    out << "}\n";
}

void write_holes_svg(const std::string &path,
                     const std::vector<Hole> &holes,
                     const SvgExportOptions &options)
{
    std::ofstream out(path, std::ios::trunc);
    require_stream(out, path);

    const Bounds b = compute_bounds(holes);
    const double w = std::max(1e-9, b.max_x - b.min_x);
    const double h = std::max(1e-9, b.max_y - b.min_y);

    const int w_px = std::max(100, options.width_px);
    const int h_px = std::max(100, options.height_px);
    const int margin_px = std::max(0, options.margin_px);

    const double usable_w = std::max(1.0, static_cast<double>(w_px - 2 * margin_px));
    const double usable_h = std::max(1.0, static_cast<double>(h_px - 2 * margin_px));
    const double scale = std::min(safe_div(usable_w, w), safe_div(usable_h, h));

    auto map_x = [&](double x) { return margin_px + (x - b.min_x) * scale; };
    auto map_y = [&](double y) {
        const double mapped = margin_px + (y - b.min_y) * scale;
        if (!options.invert_y)
        {
            return mapped;
        }
        return static_cast<double>(h_px) - mapped;
    };

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w_px << "\" height=\"" << h_px
        << "\" viewBox=\"0 0 " << w_px << " " << h_px << "\">\n";
    out << "  <rect x=\"0\" y=\"0\" width=\"" << w_px << "\" height=\"" << h_px
        << "\" fill=\"" << options.background << "\"/>\n";
    out << "  <g fill=\"none\" stroke=\"#2f3b68\" stroke-width=\"1\">\n";
    out << "    <rect x=\"" << margin_px << "\" y=\"" << margin_px << "\" width=\"" << (w_px - 2 * margin_px)
        << "\" height=\"" << (h_px - 2 * margin_px) << "\"/>\n";
    out << "  </g>\n";

    out << "  <g fill=\"" << options.fill << "\" fill-opacity=\"" << options.fill_opacity
        << "\" stroke=\"" << options.stroke << "\" stroke-opacity=\"" << options.stroke_opacity
        << "\" stroke-width=\"" << options.stroke_width << "\">\n";
    out << std::fixed << std::setprecision(3);
    for (const auto &hole : holes)
    {
        const double cx = map_x(hole.center.x());
        const double cy = map_y(hole.center.y());
        const double r = std::max(0.0, hole.radius * scale);
        out << "    <circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\"/>\n";
        if (options.draw_centers)
        {
            out << "    <circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"1.5\" fill=\"#ffcc66\"/>\n";
        }
    }
    out << "  </g>\n";
    out << "</svg>\n";
}

} // namespace photonics

