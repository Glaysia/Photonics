#ifndef PHOTONICS_SIMULATION_CONFIG_HPP
#define PHOTONICS_SIMULATION_CONFIG_HPP

#include <meep.hpp>
#include <meep/meepgeom.hpp>

#include <tuple>
#include <vector>

namespace photonics
{

// Mirror geometry aliases so users can pass/receive meep::geometric_object in
// a conventional namespace instead of the C API types.
namespace meep
{
using geometric_object = ::geometric_object;
using geometric_object_list = ::geometric_object_list;
using material_type = meep_geom::material_type;
} // namespace meep

// Parameters are expressed in units of the lattice constant a unless noted.
struct SimulationParameters
{
    double lattice_constant_um = 0.42;    // a in microns
    double resolution_px_per_a = 20.0;    // pixels per a
    double cell_size_x_in_a = 10.0;       // simulation domain size along x (in a)
    double cell_size_y_in_a = 8.0;        // simulation domain size along y (in a)
    double cell_size_z_in_a = 5.0;        // simulation domain size along z (in a)
    double pml_thickness_in_a = 1.0;      // PML thickness on each side (in a)
    double courant = 0.5;                 // CFL safety factor
};

class SimulationConfig
{
public:
    SimulationConfig();
    explicit SimulationConfig(const SimulationParameters &params);

    double lattice_constant_um() const;
    double resolution_px_per_a() const;
    double timestep() const; // dt = Courant / resolution (time unit is a / c)
    double pml_thickness_in_a() const;
    double cell_size_x_in_a() const;
    double cell_size_y_in_a() const;
    double cell_size_z_in_a() const;
    double courant() const;

    std::tuple<int, int, int> expected_grid_points() const;
    meep::grid_volume make_grid_volume() const;
    meep::boundary_region make_pml() const;
    meep::structure make_structure(const std::vector<meep::geometric_object> &geometry) const;

private:
    SimulationParameters params_;
};

} // namespace photonics

#endif // PHOTONICS_SIMULATION_CONFIG_HPP
