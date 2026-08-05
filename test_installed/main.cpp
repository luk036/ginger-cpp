#include <ginger/version.h>

#include <iostream>

auto main() -> int {
    const auto ok = (GINGER_VERSION_MAJOR >= 1);
    std::cout << "ginger installed test: version " << GINGER_VERSION << "\n";
    return ok ? 0 : 1;
}
