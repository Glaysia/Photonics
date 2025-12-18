#include "harminv_runner.hpp"

#include <meep.hpp>

#include <cmath>
#include <stdexcept>

namespace photonics
{

std::vector<HarminvMode> run_harminv(const std::vector<std::complex<double>> &samples,
                                    double dt,
                                    const HarminvRunConfig &config)
{
    if (samples.size() < 8)
    {
        return {};
    }
    if (!(dt > 0.0))
    {
        throw std::runtime_error("run_harminv: dt must be positive");
    }
    if (!(config.fmax > config.fmin))
    {
        throw std::runtime_error("run_harminv: fmax must be > fmin");
    }
    if (config.maxbands <= 0)
    {
        throw std::runtime_error("run_harminv: maxbands must be positive");
    }

    std::vector<std::complex<double>> data(samples.begin(), samples.end());
    std::vector<std::complex<double>> amps(static_cast<std::size_t>(config.maxbands));
    std::vector<double> freq_re(static_cast<std::size_t>(config.maxbands));
    std::vector<double> freq_im(static_cast<std::size_t>(config.maxbands));

    const int n = static_cast<int>(data.size());
    const int found = meep::do_harminv(data.data(),
                                       n,
                                       dt,
                                       config.fmin,
                                       config.fmax,
                                       config.maxbands,
                                       amps.data(),
                                       freq_re.data(),
                                       freq_im.data(),
                                       nullptr,
                                       1.1,
                                       config.q_threshold);

    if (found < 0)
    {
        throw std::runtime_error("run_harminv: meep::do_harminv failed");
    }

    std::vector<HarminvMode> modes;
    modes.reserve(static_cast<std::size_t>(found));
    for (int i = 0; i < found; ++i)
    {
        HarminvMode m;
        m.frequency = freq_re[static_cast<std::size_t>(i)];
        m.decay_rate = std::abs(freq_im[static_cast<std::size_t>(i)]);
        m.amplitude = std::abs(amps[static_cast<std::size_t>(i)]);
        modes.push_back(m);
    }
    return modes;
}

} // namespace photonics

