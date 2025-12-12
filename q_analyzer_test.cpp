#include "q_analyzer.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
bool approx_equal(double a, double b, double tol = 1e-6)
{
    return std::abs(a - b) <= tol;
}
} // namespace

int main()
{
    photonics::QAnalyzer analyzer;

    // Quality factor selection picks mode with largest amplitude magnitude.
    const std::vector<photonics::HarminvMode> modes = {
        {0.150, 1.0e-4, 0.2},
        {0.140, 2.0e-4, 0.9}, // dominant
        {0.130, 5.0e-4, 0.1},
    };

    const auto qres = analyzer.compute_quality(modes);
    const double expected_q = 0.140 / (2.0 * 2.0e-4);
    assert(approx_equal(qres.resonance_frequency, 0.140));
    assert(approx_equal(qres.decay_rate, 2.0e-4));
    assert(approx_equal(qres.quality_factor, expected_q));

    // Flux bookkeeping.
    const auto flux = analyzer.flux_components(1.0, 0.5, 0.25);
    assert(approx_equal(flux.total(), 1.0 + 0.5 + 0.25));

    // Mode volume check: V = U / (0.5 * eps * |E|max^2).
    const auto vol = analyzer.compute_mode_volume(3.0, 0.6, 12.0);
    const double expected_vol = 3.0 / (0.5 * 12.0 * 0.36);
    assert(approx_equal(vol.mode_volume, expected_vol));

    std::cout << "q_analyzer tests passed.\n";
    return 0;
}
