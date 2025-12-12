#include "monitor_suite.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace
{
double center_coord(const meep::volume &vol, meep::direction dir)
{
    return vol.center().in_direction(dir);
}
} // namespace

int main()
{
    photonics::MonitorConfig config;
    config.cell = {3.0, 2.0, 1.2, meep::D3};
    config.pml_thickness = 0.1;
    config.monitor_thickness = 0.02;

    photonics::MonitorSuite suite(config);
    suite.setup_default_flux_boxes();
    suite.setup_default_field_snapshot(meep::Z, 0.0);
    suite.add_harminv(meep::vec(0.0, 0.0, 0.0), 0.28, 0.08, meep::Ez);

    assert(suite.flux_monitors().size() == 6);
    assert(suite.field_monitors().size() == 1);
    assert(suite.harminv_monitors().size() == 1);

    const double expected_z = 0.5 * config.cell.sz - config.pml_thickness;
    const auto &fluxes = suite.flux_monitors();
    const auto &top = fluxes[4]; // ordering from setup_default_flux_boxes (front/back + top/bottom)
    const auto &bottom = fluxes[5];

    assert(std::fabs(center_coord(top.region, meep::Z) - expected_z) < 1e-9);
    assert(std::fabs(center_coord(bottom.region, meep::Z) + expected_z) < 1e-9);

    const auto mode = photonics::MonitorSuite::summarize_mode(0.25, 5e-4, 1.0);
    assert(std::fabs(mode.q - (0.25 / (2.0 * 5e-4))) < 1e-9);

    const auto flux = photonics::MonitorSuite::summarize_flux(1.0, 0.9, 0.1, 0.1, 0.05, 0.05);
    assert(std::fabs(flux.top - 1.0) < 1e-9);
    assert(std::fabs(flux.bottom - 0.9) < 1e-9);

    std::cout << "monitor_suite_tests passed\n";
    return 0;
}
