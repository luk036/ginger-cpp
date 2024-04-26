#include <doctest/doctest.h>
#include <ginger/greeter.h>
#include <ginger/version.h>

#include <string>

TEST_CASE("Ginger") {
  using namespace ginger;

  Ginger ginger("Tests");

  CHECK(ginger.greet(LanguageCode::EN) == "Hello, Tests!");
  CHECK(ginger.greet(LanguageCode::DE) == "Hallo Tests!");
  CHECK(ginger.greet(LanguageCode::ES) == "¡Hola Tests!");
  CHECK(ginger.greet(LanguageCode::FR) == "Bonjour Tests!");
}

TEST_CASE("Ginger version") {
  static_assert(std::string_view(GINGER_VERSION) == std::string_view("1.0"));
  CHECK(std::string(GINGER_VERSION) == std::string("1.0"));
}
