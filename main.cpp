#include <algorithm>
#include <cmath>
#include <iostream>

#include "mylib.hpp"

int main()
{
    photonics::LatticeParams params;
    params.lattice_constant = 0.42;
    params.hole_radius = 0.29 * params.lattice_constant;
    params.slab_thickness = 0.6 * params.lattice_constant;

    const auto nodes = photonics::build_triangular_lattice(5, 4, params.lattice_constant);
    const auto envelope_value = photonics::gaussian_envelope(0.0, params.lattice_constant * 0.1);

    photonics::Diagnostics diag("Initialization");
    diag.report(static_cast<double>(nodes.size()), "Lattice points");
    diag.report(round(1e6 * envelope_value) / 1e6, "Gaussian envelope (center)");

    std::cout << "Sample lattice nodes:\n";
    const auto limit = std::min(nodes.size(), static_cast<std::size_t>(5));
    for (std::size_t i = 0; i < limit; ++i)
    {
        const auto &node = nodes[i];
        std::cout << "  - node " << (i + 1) << ": (" << node.first << ", " << node.second << ")\n";
    }

    std::cout << "basic properties -> "
              << "a=" << params.lattice_constant << " μm, "
              << "R=" << params.hole_radius << " μm, "
              << "T=" << params.slab_thickness << " μm\n";

    return 0;
}
