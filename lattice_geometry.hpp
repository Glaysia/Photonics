#ifndef LATTICE_GEOMETRY_HPP
#define LATTICE_GEOMETRY_HPP

#include <cstddef>
#include <vector>

#include <meep.hpp>
#include <meep/meepgeom.hpp>
#include <ctlgeom.h>

// The installed Meep headers expose geometry types in meep_geom. Mirror the
// key aliases into the meep namespace so callers can use the more familiar
// meep::geometric_object signature.
namespace meep
{
using geometric_object = ::geometric_object;
using geometric_object_list = ::geometric_object_list;
using material_type = meep_geom::material_type;
} // namespace meep

namespace photonics
{
struct LatticeGeometryParams
{
    double lattice_constant = 1.0;            // lattice spacing in simulation units (typically a=1)
    double hole_radius = 0.29 * lattice_constant;
    double slab_thickness = 0.6 * lattice_constant;
    double delta_x1 = 0.15 * lattice_constant; // outward shift for nearest neighbors
    double delta_x2 = 0.05 * lattice_constant; // outward shift for next-nearest neighbors
    double edge_radius = hole_radius - 0.01 * lattice_constant;
    double slab_index = 3.4; // default silicon-like index
    double hole_index = 1.0; // air holes
    std::size_t columns = 11;
    std::size_t rows = 9;
};

struct Hole
{
    meep::vec center{};
    double radius = 0.0;
};

class LatticeGeometry
{
public:
    explicit LatticeGeometry(LatticeGeometryParams params = {});

    const LatticeGeometryParams &params() const noexcept;

    // Returns the positioned holes (centered on the L3 cavity) without turning
    // them into Meep geometry.
    std::vector<Hole> generate_holes() const;

    // Creates a slab block plus air holes as Meep geometric_objects.
    std::vector<meep::geometric_object> build_geometry() const;

private:
    LatticeGeometryParams params_;
};

} // namespace photonics

#endif // LATTICE_GEOMETRY_HPP
