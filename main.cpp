
#include "q_analyzer.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

int main()
{
    photonics::QAnalyzer analyzer;

    // Dummy Harminv-like modes; pick the dominant amplitude to compute Q.
    const std::vector<photonics::HarminvMode> modes = {
        {0.120, 1.5e-4, 0.2},
        {0.115, 8.0e-5, 1.0},
        {0.098, 1.0e-3, 0.05},
    };

    const auto q_res = analyzer.compute_quality(modes);
    std::cout << "Resonance frequency: " << std::fixed << std::setprecision(6) << q_res.resonance_frequency
              << " (1/a)\n";
    std::cout << "Decay rate: " << q_res.decay_rate << "\n";
    std::cout << "Quality factor: " << q_res.quality_factor << "\n";

    // Flux decomposition (top/bottom/side).
    const auto flux = analyzer.flux_components(1.2, 0.8, 0.5);
    std::cout << "Flux totals (up, down, side, total): " << flux.upward << ", " << flux.downward << ", "
              << flux.lateral << ", " << flux.total() << "\n";

    // Mode volume example using energy integral and |E|max.
    const auto volume = analyzer.compute_mode_volume(/*energy_integral=*/2.0,
                                                     /*max_field_magnitude=*/0.5,
                                                     /*relative_permittivity=*/12.0);
    std::cout << "Mode volume (a^3): " << volume.mode_volume << "\n";

    return 0;
}
