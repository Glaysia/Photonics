#ifndef PHOTONICS_SOURCE_MANAGER_HPP
#define PHOTONICS_SOURCE_MANAGER_HPP

#include <memory>
#include <vector>

#include <meep.hpp>

// Lightweight descriptor for Meep sources. Meep's public headers do not ship a
// simple POD type, so we define one in the meep namespace for convenience.
namespace meep {
struct source
{
    source()
        : where(meep::vec(0.0, 0.0, 0.0))
    {
    }

    meep::volume where;
    meep::component component = meep::Ex;
    std::shared_ptr<meep::src_time> time;
    double amplitude = 1.0;
};
} // namespace meep

namespace photonics
{
struct SourceConfig
{
    double center_frequency = 0.26;                  // in units of 1/a
    double bandwidth = 0.08;                         // gaussian fwidth in 1/a (wider to inject more energy)
    double cutoff = 3.0;                             // width multiplier for envelope taper (shorter pulse)
    double amplitude = 6.0;                          // scalar applied to all sources (stronger drive)
    std::vector<meep::component> components{meep::Ex,
                                            meep::Ey}; // TE-like defaults (in-plane E)
    std::vector<meep::vec> positions{meep::vec(0.0, 0.0, 0.0)};
};

class SourceManager
{
public:
    SourceManager();
    explicit SourceManager(SourceConfig cfg);

    const SourceConfig &config() const;
    void set_config(const SourceConfig &cfg);

    // Build sources for all configured positions/components.
    std::vector<meep::source> build_sources() const;

    // Helper to construct a point volume at the given location.
    static meep::volume point_volume(const meep::vec &position);

private:
    SourceConfig config_;
    std::shared_ptr<meep::gaussian_src_time> make_time_profile() const;
};

} // namespace photonics

#endif // PHOTONICS_SOURCE_MANAGER_HPP
