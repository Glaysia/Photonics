#include "lattice_geometry.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

namespace
{
double distance2d(double x, double y)
{
    return std::sqrt(x * x + y * y);
}

bool approx_equal(double lhs, double rhs, double tol = 1e-6)
{
    return std::abs(lhs - rhs) <= tol;
}

const photonics::Hole *find_hole_near(const std::vector<photonics::Hole> &holes,
                                      double target_x,
                                      double target_y,
                                      double tol = 1e-3)
{
    for (const auto &hole : holes)
    {
        const double dx = hole.center.x() - target_x;
        const double dy = hole.center.y() - target_y;
        if (std::abs(dx) < tol && std::abs(dy) < tol)
        {
            return &hole;
        }
    }
    return nullptr;
}
} // namespace

int main()
{
    photonics::LatticeGeometry geom{};
    const auto params = geom.params();
    const auto holes = geom.generate_holes();

    const std::size_t expected_holes = params.columns * params.rows - 3;
    assert(holes.size() == expected_holes && "hole count should drop by 3 for L3 defect");

    const auto close_to_origin = std::find_if(
        holes.begin(), holes.end(), [](const photonics::Hole &h) {
            return distance2d(h.center.x(), h.center.y()) < 0.05;
        });
    assert(close_to_origin == holes.end() && "defect region should be empty near origin");

    // For the default odd counts, the defect row is centered and has no offset.
    const double a = params.lattice_constant;
    const double expected_y = 0.0;
    const double near_left_x = -2.0 * a - params.delta_x1;
    const double near_right_x = 2.0 * a + params.delta_x1;
    const double far_left_x = -3.0 * a - params.delta_x2;
    const double far_right_x = 3.0 * a + params.delta_x2;

    const auto *left = find_hole_near(holes, near_left_x, expected_y);
    const auto *right = find_hole_near(holes, near_right_x, expected_y);
    const auto *far_left = find_hole_near(holes, far_left_x, expected_y);
    const auto *far_right = find_hole_near(holes, far_right_x, expected_y);

    assert(left && right && far_left && far_right && "adjacent holes should exist around the defect");

    assert(approx_equal(left->radius, params.edge_radius));
    assert(approx_equal(right->radius, params.edge_radius));
    assert(approx_equal(far_left->radius, params.edge_radius));
    assert(approx_equal(far_right->radius, params.edge_radius));

    const auto geometry = geom.build_geometry();
    assert(geometry.size() == holes.size() + 1 && "geometry should include slab + holes");

    std::cout << "lattice_geometry_test passed (" << holes.size() << " holes)\n";
    return 0;
}
