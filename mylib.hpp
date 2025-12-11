#ifndef MYLIB_HPP
#define MYLIB_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <mfem.hpp>

namespace photonics
{
struct LatticeParams
{
    double lattice_constant = 0.0;
    double hole_radius = 0.0;
    double slab_thickness = 0.0;
};

std::vector<std::pair<double, double>> build_triangular_lattice(std::size_t columns,
                                                                std::size_t rows,
                                                                double spacing);

double gaussian_envelope(double position, double width);

class Diagnostics
{
public:
    explicit Diagnostics(std::string name);
    void report(double value, const std::string &label) const;

private:
    std::string name_;
};

} // namespace photonics

#endif // MYLIB_HPP
