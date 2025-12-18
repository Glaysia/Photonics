#include "simulation_runner.hpp"

#include <complex>
#include <stdexcept>

namespace photonics
{
namespace
{
double unity_eps(const meep::vec &)
{
    return 1.0;
}

bool is_point_volume(const meep::volume &v)
{
    return v.diameter() <= 1e-12;
}

} // namespace

SimulationRunResult SimulationRunner::run(const SimulationConfig &config,
                                          const std::vector<meep::geometric_object> &geometry,
                                          const std::vector<meep::source> &sources,
                                          const SimulationRunConfig &run_config) const
{
    SimulationRunResult result{};
    meep::simple_material_function background(&unity_eps);
    meep::grid_volume volume = config.make_grid_volume();
    meep::boundary_region br = config.make_pml();

    result.structure = std::make_unique<meep::structure>(volume, background, br, meep::identity(), 0, config.courant());
    if (!geometry.empty())
    {
        auto owned = geometry;
        meep::geometric_object_list geom_list{};
        geom_list.num_items = static_cast<int>(owned.size());
        geom_list.items = owned.data();
        meep_geom::set_materials_from_geometry(result.structure.get(), geom_list);
    }
    result.fields = std::make_unique<meep::fields>(result.structure.get());

    for (const auto &src : sources)
    {
        if (!src.time)
        {
            throw std::runtime_error("SimulationRunner: source time profile is null");
        }
        if (!is_point_volume(src.where))
        {
            throw std::runtime_error("SimulationRunner: only point sources are supported (volume source not implemented)");
        }

        result.fields->add_point_source(src.component,
                                        *src.time,
                                        src.where.center(),
                                        std::complex<double>(src.amplitude, 0.0));
    }

    result.last_source_time = result.fields->last_source_time();
    result.stop_time = result.last_source_time + run_config.until_after_sources;
    result.steps = 0;

    while (result.fields->time() + 0.5 * result.fields->dt < result.stop_time)
    {
        result.fields->step();
        ++result.steps;
    }

    return result;
}

} // namespace photonics
