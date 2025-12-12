#include "q_analyzer.hpp"

#include <algorithm>
#include <cmath>

namespace photonics
{

QualityResult QAnalyzer::compute_quality(const std::vector<HarminvMode> &modes) const
{
    QualityResult result{};
    if (modes.empty())
    {
        return result;
    }

    // Pick mode with largest amplitude magnitude as the dominant resonance.
    const auto dominant = std::max_element(
        modes.begin(), modes.end(), [](const HarminvMode &a, const HarminvMode &b) {
            return std::abs(a.amplitude) < std::abs(b.amplitude);
        });

    const double decay = std::abs(dominant->decay_rate);
    const double freq = dominant->frequency;

    result.resonance_frequency = freq;
    result.decay_rate = decay;
    result.quality_factor = (decay > kEpsilon) ? freq / (2.0 * decay) : std::numeric_limits<double>::infinity();
    return result;
}

FluxBreakdown QAnalyzer::flux_components(double upward, double downward, double lateral) const
{
    FluxBreakdown breakdown{upward, downward, lateral};
    return breakdown;
}

VolumeResult QAnalyzer::compute_mode_volume(double energy_integral,
                                            double max_field_magnitude,
                                            double relative_permittivity) const
{
    VolumeResult vol{};
    vol.energy_integral = energy_integral;
    vol.max_field_magnitude = max_field_magnitude;
    vol.relative_permittivity = relative_permittivity;

    const double denom = 0.5 * relative_permittivity * max_field_magnitude * max_field_magnitude;
    vol.mode_volume = safe_divide(energy_integral, denom);
    return vol;
}

double QAnalyzer::safe_divide(double numerator, double denominator)
{
    if (std::abs(denominator) < kEpsilon)
    {
        return std::numeric_limits<double>::infinity();
    }
    return numerator / denominator;
}

} // namespace photonics
