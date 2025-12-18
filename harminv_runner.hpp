#ifndef PHOTONICS_HARMINV_RUNNER_HPP
#define PHOTONICS_HARMINV_RUNNER_HPP

#include "q_analyzer.hpp"

#include <complex>
#include <vector>

namespace photonics
{

struct HarminvRunConfig
{
    double fmin = 0.0;
    double fmax = 0.0;
    int maxbands = 50;
    double q_threshold = 50.0;
};

// Runs Meep's harmonic inversion (Harminv) on a uniformly sampled complex time-series.
// - `dt` is the sample spacing in Meep time units.
// - `fmin`/`fmax` are in Meep frequency units (same units used for sources).
std::vector<HarminvMode> run_harminv(const std::vector<std::complex<double>> &samples,
                                    double dt,
                                    const HarminvRunConfig &config);

} // namespace photonics

#endif // PHOTONICS_HARMINV_RUNNER_HPP

