
#include "lattice_geometry.hpp"

#include <iomanip>
#include <iostream>

int main()
{
    photonics::LatticeGeometry geom{};
    const auto holes = geom.generate_holes();
    const auto geometry = geom.build_geometry();

    std::cout << "Lattice geometry built with "
              << holes.size() << " holes and "
              << geometry.size() << " geometric objects (including slab)"
              << "\n";

    const auto params = geom.params();
    std::cout << "Defaults: a=" << params.lattice_constant << " μm, "
              << "r=" << params.hole_radius << " μm, "
              << "t=" << params.slab_thickness << " μm, "
              << "Δx1=" << params.delta_x1 << " μm, "
              << "Δx2=" << params.delta_x2 << " μm"
              << "\n";

    std::cout << "First five hole centers (x, y, r):\n";
    const std::size_t limit = std::min<std::size_t>(5, holes.size());
    for (std::size_t i = 0; i < limit; ++i)
    {
        const auto &h = holes[i];
        std::cout << "  [" << (i + 1) << "] "
                  << std::fixed << std::setprecision(6)
                  << h.center.x() << ", " << h.center.y()
                  << " (r=" << h.radius << ")"
                  << "\n";
    }

    return 0;
}
