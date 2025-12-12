#include "lattice_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace photonics
{
namespace
{
constexpr double sqrt3_over_2()
{
    return 0.86602540378443864676372317075294; // std::sqrt(3.0) / 2.0;
}

meep::vec subtract(const meep::vec &lhs, const meep::vec &rhs)
{
    return meep::vec{lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z()};
}

vector3 to_vector3(const meep::vec &v)
{
    return meep_geom::make_vector3(v.x(), v.y(), v.z());
}
} // namespace

LatticeGeometry::LatticeGeometry(LatticeGeometryParams params)
    : params_(params)
{
}

const LatticeGeometryParams &LatticeGeometry::params() const noexcept
{
    return params_;
}

std::vector<Hole> LatticeGeometry::generate_holes() const
{
    const double spacing = params_.lattice_constant;
    const double offset_magnitude = spacing * 0.5;
    const std::size_t defect_row = params_.rows / 2;
    const std::size_t defect_center_col = params_.columns / 2;

    std::vector<Hole> holes;
    holes.reserve(params_.rows * params_.columns);

    std::vector<meep::vec> removed;
    removed.reserve(3);

    for (std::size_t row = 0; row < params_.rows; ++row)
    {
        const double y = static_cast<double>(row) * spacing * sqrt3_over_2();
        const bool odd_row = (row % 2 != 0);
        const double row_offset = odd_row ? offset_magnitude : 0.0;

        for (std::size_t col = 0; col < params_.columns; ++col)
        {
            double x = static_cast<double>(col) * spacing + row_offset;
            double radius = params_.hole_radius;

            const bool in_defect_band = (row == defect_row &&
                                         col >= defect_center_col - 1 &&
                                         col <= defect_center_col + 1);
            if (in_defect_band)
            {
                removed.emplace_back(meep::vec{x, y, 0.0});
                continue;
            }

            const bool in_defect_row = (row == defect_row);
            if (in_defect_row)
            {
                if (col == defect_center_col - 2)
                {
                    x -= params_.delta_x1;
                    radius = params_.edge_radius;
                }
                else if (col == defect_center_col + 2)
                {
                    x += params_.delta_x1;
                    radius = params_.edge_radius;
                }
                else if (col == defect_center_col - 3)
                {
                    x -= params_.delta_x2;
                    radius = params_.edge_radius;
                }
                else if (col == defect_center_col + 3)
                {
                    x += params_.delta_x2;
                    radius = params_.edge_radius;
                }
            }

            holes.push_back(Hole{meep::vec{x, y, 0.0}, radius});
        }
    }

    meep::vec cavity_center{0.0, 0.0, 0.0};
    if (!removed.empty())
    {
        const double sum_x = std::accumulate(removed.begin(), removed.end(), 0.0,
                                             [](double acc, const meep::vec &v) { return acc + v.x(); });
        const double sum_y = std::accumulate(removed.begin(), removed.end(), 0.0,
                                             [](double acc, const meep::vec &v) { return acc + v.y(); });
        const double inv_count = 1.0 / static_cast<double>(removed.size());
        cavity_center = meep::vec{sum_x * inv_count, sum_y * inv_count, 0.0};
    }

    for (auto &hole : holes)
    {
        hole.center = subtract(hole.center, cavity_center);
    }

    return holes;
}

std::vector<meep::geometric_object> LatticeGeometry::build_geometry() const
{
    const auto holes = generate_holes();
    std::vector<meep::geometric_object> geometry;
    geometry.reserve(holes.size() + 1);

    // Estimate a slab that comfortably encompasses all holes plus a padding.
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const auto &hole : holes)
    {
        min_x = std::min(min_x, hole.center.x());
        max_x = std::max(max_x, hole.center.x());
        min_y = std::min(min_y, hole.center.y());
        max_y = std::max(max_y, hole.center.y());
    }

    const double padding = params_.lattice_constant;
    const double slab_x = (max_x - min_x) + 2.0 * padding;
    const double slab_y = (max_y - min_y) + 2.0 * padding;

    const meep::material_type slab_material =
        meep_geom::make_dielectric(params_.slab_index * params_.slab_index);
    const meep::material_type air_material =
        meep_geom::make_dielectric(params_.hole_index * params_.hole_index);

    geometry.push_back(::make_block(slab_material,
                                    meep_geom::make_vector3(0.0, 0.0, 0.0),
                                    meep_geom::make_vector3(1.0, 0.0, 0.0),
                                    meep_geom::make_vector3(0.0, 1.0, 0.0),
                                    meep_geom::make_vector3(0.0, 0.0, 1.0),
                                    meep_geom::make_vector3(slab_x, slab_y, params_.slab_thickness)));

    for (const auto &hole : holes)
    {
        geometry.push_back(::make_cylinder(
            air_material, to_vector3(hole.center), hole.radius, params_.slab_thickness,
            meep_geom::make_vector3(0.0, 0.0, 1.0)));
    }

    return geometry;
}

} // namespace photonics
