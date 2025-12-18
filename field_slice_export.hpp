#ifndef PHOTONICS_FIELD_SLICE_EXPORT_HPP
#define PHOTONICS_FIELD_SLICE_EXPORT_HPP

#include <meep.hpp>

#include <string>

namespace photonics
{

struct FieldSliceExportConfig
{
    meep::component component = meep::Ey;
    double z = 0.0;

    // Output grid size (samples along x/y).
    int nx = 0;
    int ny = 0;

    // Physical extents in simulation coordinates (typically in units of a).
    double sx = 0.0;
    double sy = 0.0;

    // When true, exports the imaginary part instead of the real part.
    bool imag_part = false;
};

// Samples a field component on a z=const plane and writes a 2D CSV matrix.
// Format:
//   - Header comment line with metadata
//   - ny rows, each containing nx comma-separated values
void write_field_slice_csv(const std::string &path, meep::fields &fields, const FieldSliceExportConfig &cfg);

} // namespace photonics

#endif // PHOTONICS_FIELD_SLICE_EXPORT_HPP

