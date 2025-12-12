#ifndef PHOTONICS_Q_ANALYZER_HPP
#define PHOTONICS_Q_ANALYZER_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace photonics
{
// Simple DTO representing a Harminv-like mode entry.
struct HarminvMode
{
    double frequency = 0.0;   // Resonant frequency (in units of 1/a).
    double decay_rate = 0.0;  // Exponential decay rate (imag part), positive for loss.
    double amplitude = 0.0;   // Optional weighting for selecting the dominant mode.
};

// Resulting resonance metrics.
struct QualityResult
{
    double resonance_frequency = 0.0;
    double decay_rate = 0.0;
    double quality_factor = 0.0;
};

// Flux bookkeeping for top/bottom/side channels.
struct FluxBreakdown
{
    double upward = 0.0;
    double downward = 0.0;
    double lateral = 0.0;
    double total() const { return upward + downward + lateral; }
};

// Mode volume calculation result.
struct VolumeResult
{
    double mode_volume = std::numeric_limits<double>::infinity();
    double energy_integral = 0.0;
    double max_field_magnitude = 0.0;
    double relative_permittivity = 1.0;
};

class QAnalyzer
{
public:
    QAnalyzer() = default;

    QualityResult compute_quality(const std::vector<HarminvMode> &modes) const;
    FluxBreakdown flux_components(double upward, double downward, double lateral) const;
    VolumeResult compute_mode_volume(double energy_integral,
                                     double max_field_magnitude,
                                     double relative_permittivity = 1.0) const;

private:
    static constexpr double kEpsilon = 1e-12;
    static double safe_divide(double numerator, double denominator);
};

} // namespace photonics

#endif // PHOTONICS_Q_ANALYZER_HPP
