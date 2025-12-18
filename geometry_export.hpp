#ifndef PHOTONICS_GEOMETRY_EXPORT_HPP
#define PHOTONICS_GEOMETRY_EXPORT_HPP

#include "lattice_geometry.hpp"

#include <string>
#include <vector>

namespace photonics
{

struct GeometryMetadata
{
    LatticeGeometryParams params{};
    std::vector<Hole> holes{};
    std::string units = "um";
};

struct SvgExportOptions
{
    int width_px = 1200;
    int height_px = 900;
    int margin_px = 60;
    bool draw_centers = false;
    bool invert_y = true;
    std::string background = "#0b1020";
    std::string stroke = "#dfe7ff";
    std::string fill = "#dfe7ff";
    double fill_opacity = 0.10;
    double stroke_opacity = 0.85;
    double stroke_width = 1.0;
};

void write_holes_csv(const std::string &path, const std::vector<Hole> &holes);
void write_geometry_json(const std::string &path, const GeometryMetadata &metadata);
void write_holes_svg(const std::string &path,
                     const std::vector<Hole> &holes,
                     const SvgExportOptions &options = {});

} // namespace photonics

#endif // PHOTONICS_GEOMETRY_EXPORT_HPP

