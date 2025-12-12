#include "monitor_suite.hpp"

#include <algorithm>
#include <cmath>

namespace photonics
{

namespace
{
double half_extent(double full_size)
{
    return 0.5 * full_size;
}
} // namespace

MonitorSuite::MonitorSuite(const MonitorConfig &config)
    : config_(config)
{
}

const MonitorConfig &MonitorSuite::config() const
{
    return config_;
}

HarminvMonitor MonitorSuite::add_harminv(const meep::vec &point,
                                         double center_frequency,
                                         double bandwidth,
                                         meep::component component)
{
    HarminvMonitor monitor{point, center_frequency, bandwidth, component};
    harminv_monitors_.push_back(monitor);
    return monitor;
}

FluxMonitor MonitorSuite::add_flux_box(const std::string &label,
                                       meep::direction normal,
                                       double coordinate,
                                       double thickness)
{
    FluxMonitor monitor{label, build_face(normal, coordinate, thickness)};
    flux_monitors_.push_back(monitor);
    return monitor;
}

FieldSnapshot MonitorSuite::add_field_snapshot(const std::string &label,
                                               meep::direction normal,
                                               double coordinate,
                                               double thickness,
                                               meep::component component)
{
    FieldSnapshot snapshot{label, build_face(normal, coordinate, thickness), component};
    field_monitors_.push_back(snapshot);
    return snapshot;
}

void MonitorSuite::setup_default_flux_boxes()
{
    flux_monitors_.clear();

    const double thickness = std::max(1e-4, config_.monitor_thickness);
    const double sx = config_.cell.sx;
    const double sy = config_.cell.sy;
    const double sz = config_.cell.sz;

    const double x_face = half_extent(sx) - config_.pml_thickness;
    const double y_face = half_extent(sy) - config_.pml_thickness;
    const double z_face = half_extent(sz) - config_.pml_thickness;

    add_flux_box("flux-left", meep::X, -x_face, thickness);
    add_flux_box("flux-right", meep::X, x_face, thickness);
    add_flux_box("flux-front", meep::Y, y_face, thickness);
    add_flux_box("flux-back", meep::Y, -y_face, thickness);

    if (config_.cell.dim == meep::D3)
    {
        add_flux_box("flux-top", meep::Z, z_face, thickness);
        add_flux_box("flux-bottom", meep::Z, -z_face, thickness);
    }
    else
    {
        // 2D case still exposes top/bottom as Y-facing slabs for compatibility.
        add_flux_box("flux-top", meep::Y, y_face, thickness);
        add_flux_box("flux-bottom", meep::Y, -y_face, thickness);
    }
}

void MonitorSuite::setup_default_field_snapshot(meep::direction normal, double coordinate)
{
    const double thickness = std::max(1e-4, config_.monitor_thickness);
    add_field_snapshot("field-snapshot", normal, coordinate, thickness);
}

const std::vector<HarminvMonitor> &MonitorSuite::harminv_monitors() const
{
    return harminv_monitors_;
}

const std::vector<FluxMonitor> &MonitorSuite::flux_monitors() const
{
    return flux_monitors_;
}

const std::vector<FieldSnapshot> &MonitorSuite::field_monitors() const
{
    return field_monitors_;
}

HarminvResult MonitorSuite::summarize_mode(double frequency, double decay_rate, double amplitude)
{
    HarminvResult result;
    result.frequency = frequency;
    result.decay_rate = decay_rate;
    result.amplitude = amplitude;
    result.q = decay_rate > 0.0 ? frequency / (2.0 * decay_rate) : 0.0;
    return result;
}

FluxResult MonitorSuite::summarize_flux(double top,
                                        double bottom,
                                        double left,
                                        double right,
                                        double front,
                                        double back)
{
    FluxResult result;
    result.top = top;
    result.bottom = bottom;
    result.left = left;
    result.right = right;
    result.front = front;
    result.back = back;
    return result;
}

meep::volume MonitorSuite::build_face(meep::direction normal, double coordinate, double thickness) const
{
    const double sx = config_.cell.sx;
    const double sy = config_.cell.sy;
    const double sz = config_.cell.sz;
    const double half_x = half_extent(sx);
    const double half_y = half_extent(sy);
    const double half_z = config_.cell.dim == meep::D3 ? half_extent(sz) : 0.0;
    const double half_t = half_extent(thickness);

    if (config_.cell.dim == meep::D2)
    {
        if (normal == meep::X)
        {
            return meep::volume(meep::vec(coordinate - half_t, -half_y),
                                meep::vec(coordinate + half_t, half_y));
        }
        return meep::volume(meep::vec(-half_x, coordinate - half_t),
                            meep::vec(half_x, coordinate + half_t));
    }

    switch (normal)
    {
    case meep::X:
        return meep::volume(meep::vec(coordinate - half_t, -half_y, -half_z),
                            meep::vec(coordinate + half_t, half_y, half_z));
    case meep::Y:
        return meep::volume(meep::vec(-half_x, coordinate - half_t, -half_z),
                            meep::vec(half_x, coordinate + half_t, half_z));
    case meep::Z:
    default:
        return meep::volume(meep::vec(-half_x, -half_y, coordinate - half_t),
                            meep::vec(half_x, half_y, coordinate + half_t));
    }
}

} // namespace photonics
