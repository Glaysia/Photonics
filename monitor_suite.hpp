#ifndef MONITOR_SUITE_HPP
#define MONITOR_SUITE_HPP

#include <meep.hpp>

#include <string>
#include <vector>

namespace photonics
{

struct CellDimensions
{
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    meep::ndim dim = meep::D3;
};

struct MonitorConfig
{
    CellDimensions cell;
    double pml_thickness = 0.0;
    double monitor_thickness = 0.02;
};

struct HarminvResult
{
    double frequency = 0.0;
    double decay_rate = 0.0;
    double q = 0.0;
    double amplitude = 0.0;
};

struct FluxResult
{
    double top = 0.0;
    double bottom = 0.0;
    double left = 0.0;
    double right = 0.0;
    double front = 0.0;
    double back = 0.0;
};

struct FieldSnapshot
{
    std::string label;
    meep::volume region;
    meep::component component = meep::Ez;
};

struct FluxMonitor
{
    std::string label;
    meep::volume region;
};

struct HarminvMonitor
{
    meep::vec location;
    double center_frequency = 0.0;
    double bandwidth = 0.0;
    meep::component component = meep::Ez;
};

class MonitorSuite
{
public:
    explicit MonitorSuite(const MonitorConfig &config);

    const MonitorConfig &config() const;

    HarminvMonitor add_harminv(const meep::vec &point,
                               double center_frequency,
                               double bandwidth,
                               meep::component component = meep::Ez);

    FluxMonitor add_flux_box(const std::string &label,
                             meep::direction normal,
                             double coordinate,
                             double thickness);

    FieldSnapshot add_field_snapshot(const std::string &label,
                                     meep::direction normal,
                                     double coordinate,
                                     double thickness,
                                     meep::component component = meep::Ez);

    void setup_default_flux_boxes();
    void setup_default_field_snapshot(meep::direction normal = meep::Z, double coordinate = 0.0);

    const std::vector<HarminvMonitor> &harminv_monitors() const;
    const std::vector<FluxMonitor> &flux_monitors() const;
    const std::vector<FieldSnapshot> &field_monitors() const;

    static HarminvResult summarize_mode(double frequency, double decay_rate, double amplitude);
    static FluxResult summarize_flux(double top,
                                     double bottom,
                                     double left,
                                     double right,
                                     double front,
                                     double back);

private:
    meep::volume build_face(meep::direction normal, double coordinate, double thickness) const;

    MonitorConfig config_;
    std::vector<HarminvMonitor> harminv_monitors_;
    std::vector<FluxMonitor> flux_monitors_;
    std::vector<FieldSnapshot> field_monitors_;
};

} // namespace photonics

#endif // MONITOR_SUITE_HPP
