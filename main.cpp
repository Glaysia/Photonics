
#include "monitor_suite.hpp"

#include <iomanip>
#include <iostream>

int main()
{
    const double a = 0.42;          // lattice constant (μm)
    const double slab_thickness = 0.6 * a;
    const double pml = 0.2 * a;

    photonics::MonitorConfig config;
    config.cell = {5.0 * a, 5.0 * a, slab_thickness + 0.4 * a, meep::D3};
    config.pml_thickness = pml;
    config.monitor_thickness = 0.05 * a;

    photonics::MonitorSuite suite(config);
    suite.setup_default_flux_boxes();
    suite.setup_default_field_snapshot(meep::Z, 0.0);
    suite.add_harminv(meep::vec(0.0, 0.0, 0.0), 0.25, 0.05, meep::Ez);

    std::cout << "Monitor suite initialized with cell sizes (μm): "
              << config.cell.sx << " x " << config.cell.sy << " x " << config.cell.sz << "\n";
    std::cout << "Flux monitor count: " << suite.flux_monitors().size() << "\n";
    std::cout << "Field snapshots: " << suite.field_monitors().size()
              << ", Harminv probes: " << suite.harminv_monitors().size() << "\n";

    for (const auto &flux : suite.flux_monitors())
    {
        const auto center = flux.region.center();
        std::cout << "  - " << flux.label << " center (μm) -> "
                  << "x=" << std::fixed << std::setprecision(4) << center.in_direction(meep::X) << ", "
                  << "y=" << center.in_direction(meep::Y) << ", "
                  << "z=" << center.in_direction(meep::Z) << "\n";
    }

    const auto mode_summary = photonics::MonitorSuite::summarize_mode(0.25, 1e-3, 1.0);
    const auto flux_summary = photonics::MonitorSuite::summarize_flux(1.0, 0.9, 0.1, 0.1, 0.05, 0.05);

    std::cout << "Mode summary -> freq=" << mode_summary.frequency << " decay=" << mode_summary.decay_rate
              << " Q=" << mode_summary.q << "\n";
    std::cout << "Flux summary (top/bottom/lat): "
              << flux_summary.top << " / " << flux_summary.bottom << " / "
              << (flux_summary.left + flux_summary.right + flux_summary.front + flux_summary.back) << "\n";

    return 0;
}
