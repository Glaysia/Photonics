#include "source_manager.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace
{
double constant_eps(const meep::vec &)
{
    return 12.0;
}

bool approximately_equal(double a, double b, double tol = 1e-6)
{
    const double diff = std::abs(a - b);
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return diff < tol * scale;
}
} // namespace

void test_default_configuration()
{
    std::cout << "Running test_default_configuration...\n";
    photonics::SourceManager manager;
    auto sources = manager.build_sources();

    assert(sources.size() == 2); // Ex/Ey at the origin

    const auto *time = dynamic_cast<meep::gaussian_src_time *>(sources.front().time.get());
    assert(time != nullptr);
    std::cout << "  observed center_frequency=" << time->frequency().real()
              << ", configured=" << manager.config().center_frequency << "\n";
    const double observed_fwidth = time->get_fwidth();
    std::cout << "  observed fwidth=" << observed_fwidth << ", configured=" << manager.config().bandwidth << "\n";
    assert(approximately_equal(time->frequency().real(), manager.config().center_frequency));
    assert(observed_fwidth > 0.0);
    assert(observed_fwidth > manager.config().bandwidth * 0.1);
    assert(observed_fwidth < manager.config().bandwidth * 10.0);

    const auto center = sources.front().where.center();
    assert(approximately_equal(center.x(), 0.0));
    assert(approximately_equal(center.y(), 0.0));
    assert(approximately_equal(center.z(), 0.0));
}

void test_custom_configuration()
{
    std::cout << "Running test_custom_configuration...\n";
    photonics::SourceConfig cfg;
    cfg.center_frequency = 0.18;
    cfg.bandwidth = 0.07;
    cfg.cutoff = 7.0;
    cfg.amplitude = 0.75;
    cfg.components = {meep::Ey};
    cfg.positions = {meep::vec(0.5, -0.1, 0.0), meep::vec(-0.25, 0.35, 0.0)};

    photonics::SourceManager manager(cfg);
    auto sources = manager.build_sources();

    assert(sources.size() == cfg.components.size() * cfg.positions.size());

    const auto *time = dynamic_cast<meep::gaussian_src_time *>(sources.front().time.get());
    assert(time != nullptr);
    const double observed_bandwidth = time->get_fwidth();
    std::cout << "  custom observed fwidth=" << observed_bandwidth << ", configured=" << cfg.bandwidth << "\n";
    assert(approximately_equal(time->frequency().real(), cfg.center_frequency));
    assert(observed_bandwidth > cfg.bandwidth * 0.1);
    assert(observed_bandwidth < cfg.bandwidth * 10.0);

    for (std::size_t i = 0; i < sources.size(); ++i)
    {
        const auto center = sources[i].where.center();
        assert(approximately_equal(center.x(), cfg.positions[i].x()));
        assert(approximately_equal(center.y(), cfg.positions[i].y()));
        assert(approximately_equal(center.z(), cfg.positions[i].z()));
        assert(approximately_equal(sources[i].amplitude, cfg.amplitude));
        assert(sources[i].component == cfg.components.front());
    }
}

void test_structure_wiring()
{
    std::cout << "Running test_structure_wiring...\n";
    meep::grid_volume cell = meep::vol3d(1.0, 1.0, 0.6, 10.0);
    meep::structure structure(cell, constant_eps);

    photonics::SourceManager manager;
    auto sources = manager.build_sources();

    assert(cell.nx() > 0);
    assert(cell.nz() > 0);
    assert(!sources.empty());
    (void)structure; // ensure the object is instantiated and linked
}

int main()
{
    std::cout << "Starting SourceManager tests...\n";
    int argc = 1;
    char arg0[] = "source_manager_tests";
    char *argv_raw[] = {arg0, nullptr};
    char **argv = argv_raw;
    meep::initialize mpi(argc, argv);

    test_default_configuration();
    test_custom_configuration();
    test_structure_wiring();

    std::cout << "SourceManager tests passed.\n";
    return 0;
}
