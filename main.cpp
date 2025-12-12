#include "source_manager.hpp"

#include <iostream>
#include <string>

namespace
{
double uniform_slab_eps(const meep::vec &)
{
    return 12.0;
}
} // namespace

int main()
{
    int argc = 1;
    char arg0[] = "photonics_source_main";
    char *argv_raw[] = {arg0, nullptr};
    char **argv = argv_raw;
    meep::initialize mpi(argc, argv);

    meep::grid_volume cell = meep::vol3d(2.0, 2.0, 0.6, 12.0); // 12 px/a resolution
    meep::structure structure(cell, uniform_slab_eps);

    photonics::SourceManager manager;
    auto sources = manager.build_sources();

    if (sources.size() != 2)
    {
        std::cerr << "Unexpected default source count: " << sources.size() << "\n";
        return 1;
    }

    std::cout << "Structure grid (nx, ny, nz): " << cell.nx() << ", " << cell.ny() << ", " << cell.nz() << "\n";
    std::cout << "Generated " << sources.size() << " sources\n";

    const auto *time = dynamic_cast<meep::gaussian_src_time *>(sources.front().time.get());
    if (time != nullptr)
    {
        std::cout << "Center frequency: " << time->frequency().real()
                  << " | bandwidth: " << time->get_fwidth() << "\n";
    }

    std::size_t index = 0;
    for (const auto &src : sources)
    {
        const auto center = src.where.center();
        std::cout << "  [" << index++ << "] component=" << src.component << " pos=(" << center.x() << ", "
                  << center.y() << ", " << center.z() << "), amplitude=" << src.amplitude << "\n";
    }

    return 0;
}
