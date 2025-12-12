#include "simulation_config.hpp"

#include <cmath>

namespace
{
double unity_eps(const meep::vec &)
{
    return 1.0; // background material (air)
}
} // namespace

namespace photonics
{

SimulationConfig::SimulationConfig() = default;

SimulationConfig::SimulationConfig(const SimulationParameters &params)
    : params_(params)
{
}

double SimulationConfig::lattice_constant_um() const
{
    return params_.lattice_constant_um;
}

double SimulationConfig::resolution_px_per_a() const
{
    return params_.resolution_px_per_a;
}

double SimulationConfig::timestep() const
{
    return params_.courant / params_.resolution_px_per_a;
}

double SimulationConfig::pml_thickness_in_a() const
{
    return params_.pml_thickness_in_a;
}

double SimulationConfig::cell_size_x_in_a() const
{
    return params_.cell_size_x_in_a;
}

double SimulationConfig::cell_size_y_in_a() const
{
    return params_.cell_size_y_in_a;
}

double SimulationConfig::cell_size_z_in_a() const
{
    return params_.cell_size_z_in_a;
}

double SimulationConfig::courant() const
{
    return params_.courant;
}

std::tuple<int, int, int> SimulationConfig::expected_grid_points() const
{
    const int nx = static_cast<int>(std::lround(params_.cell_size_x_in_a * params_.resolution_px_per_a));
    const int ny = static_cast<int>(std::lround(params_.cell_size_y_in_a * params_.resolution_px_per_a));
    const int nz = static_cast<int>(std::lround(params_.cell_size_z_in_a * params_.resolution_px_per_a));
    return {nx, ny, nz};
}

meep::grid_volume SimulationConfig::make_grid_volume() const
{
    return meep::vol3d(params_.cell_size_x_in_a,
                       params_.cell_size_y_in_a,
                       params_.cell_size_z_in_a,
                       params_.resolution_px_per_a);
}

meep::boundary_region SimulationConfig::make_pml() const
{
    return meep::pml(params_.pml_thickness_in_a);
}

meep::structure SimulationConfig::make_structure(const std::vector<meep::geometric_object> &geometry) const
{
    meep::simple_material_function background(&unity_eps);
    meep::grid_volume volume = make_grid_volume();
    meep::boundary_region br = make_pml();

    meep::structure structure(volume, background, br, meep::identity(), 0, params_.courant);

    if (!geometry.empty())
    {
        auto owned = geometry;
        meep::geometric_object_list geom_list{};
        geom_list.num_items = static_cast<int>(owned.size());
        geom_list.items = owned.data();
        meep_geom::set_materials_from_geometry(&structure, geom_list);
    }

    return structure;
}

} // namespace photonics
