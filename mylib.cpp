#include "mylib.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

namespace photonics
{

Diagnostics::Diagnostics(std::string name)
    : name_(std::move(name))
{
}

void Diagnostics::report(double value, const std::string &label) const
{
    std::cout << "[" << name_ << "] " << label << ": " << std::fixed << std::setprecision(6) << value << "\n";
}

std::vector<std::pair<double, double>> build_triangular_lattice(std::size_t columns,
                                                                std::size_t rows,
                                                                double spacing)
{
    std::vector<std::pair<double, double>> nodes;
    nodes.reserve(columns * rows);

    constexpr double constexpr_sqrt_3_over_2 = 0.86602540378443864676372317075294;

    for (std::size_t row = 0; row < rows; ++row)
    {
        const double y = static_cast<double>(row) * spacing * constexpr_sqrt_3_over_2;
        const double offset = (row % 2 == 0 ? 0.0 : spacing / 2.0);

        for (std::size_t column = 0; column < columns; ++column)
        {
            const double x = static_cast<double>(column) * spacing + offset;
            nodes.emplace_back(x, y);
        }
    }

    return nodes;
}

double gaussian_envelope(double position, double width)
{
    if (width <= 0.0)
    {
        return position == 0.0 ? 1.0 : 0.0;
    }

    const double normalized = position / width;
    return std::exp(-0.5 * normalized * normalized);
}

} // namespace photonics
