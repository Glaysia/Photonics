#include "simulation_config.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <tuple>

int main()
{
    photonics::SimulationParameters params; // defaults anchored to a = 0.42 µm
    photonics::SimulationConfig config(params);

    const meep::grid_volume gv = config.make_grid_volume();
    const auto [expected_nx, expected_ny, expected_nz] = config.expected_grid_points();
    assert(gv.nx() == expected_nx);
    assert(gv.ny() == expected_ny);
    assert(gv.nz() == expected_nz);

    const auto pml_region = config.make_pml();
    assert(pml_region.check_ok(gv));

    const meep::structure structure = config.make_structure({});
    assert(std::abs(structure.dt - config.timestep()) < 1e-12);
    assert(std::abs(structure.Courant - config.courant()) < 1e-12);

    const double ax_um = config.cell_size_x_in_a() * config.lattice_constant_um();
    const double ay_um = config.cell_size_y_in_a() * config.lattice_constant_um();
    const double az_um = config.cell_size_z_in_a() * config.lattice_constant_um();

    std::cout << "SimulationConfig defaults (units: a; physical units shown in µm):\n";
    std::cout << "  - lattice constant a: " << config.lattice_constant_um() << " µm\n";
    std::cout << "  - resolution: " << config.resolution_px_per_a() << " px/a\n";
    std::cout << "  - cell (a): " << config.cell_size_x_in_a() << " x " << config.cell_size_y_in_a()
              << " x " << config.cell_size_z_in_a() << "\n";
    std::cout << "  - cell (µm): " << ax_um << " x " << ay_um << " x " << az_um << "\n";
    std::cout << "  - PML thickness: " << config.pml_thickness_in_a() << " a per side\n";
    std::cout << "  - timestep dt: " << config.timestep() << " (a/c units)\n";

    std::cout << "grid_volume points: " << gv.nx() << " x " << gv.ny() << " x " << gv.nz() << "\n";

    const auto freqs = meep::linspace(0.10, 0.14, 5);
    std::cout << "linspace sample:";
    for (const auto f : freqs)
    {
        std::cout << " " << f;
    }
    std::cout << "\n";

    std::cout << "All assertions passed for SimulationConfig.\n";
    return 0;
}
