#include "source_manager.hpp"

#include <algorithm>
#include <cassert>
#include <utility>

namespace photonics
{

SourceManager::SourceManager() = default;

SourceManager::SourceManager(SourceConfig cfg)
    : config_(std::move(cfg))
{
    if (config_.components.empty())
    {
        config_.components = {meep::Ex, meep::Ey};
    }
    if (config_.positions.empty())
    {
        config_.positions = {meep::vec(0.0, 0.0, 0.0)};
    }
}

const SourceConfig &SourceManager::config() const
{
    return config_;
}

void SourceManager::set_config(const SourceConfig &cfg)
{
    config_ = cfg;
}

std::shared_ptr<meep::gaussian_src_time> SourceManager::make_time_profile() const
{
    const double fwidth = std::max(config_.bandwidth, 1e-9);
    auto profile = std::make_shared<meep::gaussian_src_time>(config_.center_frequency,
                                                             fwidth,
                                                             config_.cutoff);
    profile->is_integrated = true; // TE-like electric source by default
    return profile;
}

meep::volume SourceManager::point_volume(const meep::vec &position)
{
    return meep::volume(position);
}

std::vector<meep::source> SourceManager::build_sources() const
{
    std::vector<meep::source> sources;
    const auto profile = make_time_profile();

    const std::vector<meep::component> components =
        config_.components.empty() ? std::vector<meep::component>{meep::Ex, meep::Ey}
                                   : config_.components;
    const std::vector<meep::vec> positions =
        config_.positions.empty() ? std::vector<meep::vec>{meep::vec(0.0, 0.0, 0.0)}
                                  : config_.positions;

    for (const auto &pos : positions)
    {
        const meep::volume where = point_volume(pos);
        for (const auto comp : components)
        {
            meep::source s;
            s.where = where;
            s.component = comp;
            s.time = profile;
            s.amplitude = config_.amplitude;
            sources.push_back(std::move(s));
        }
    }

    return sources;
}

} // namespace photonics
