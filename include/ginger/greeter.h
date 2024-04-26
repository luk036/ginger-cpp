#pragma once

#include <string>

namespace ginger {

  /**  Language codes to be used with the Ginger class */
  enum class LanguageCode { EN, DE, ES, FR };

  /**
   * @brief A class for saying hello in multiple languages
   */
  class Ginger {
    std::string name;

  public:
    /**
     * @brief Creates a new ginger
     * @param name the name to greet
     */
    Ginger(std::string name);

    /**
     * @brief Creates a localized string containing the greeting
     * @param lang the language to greet in
     * @return a string containing the greeting
     */
    std::string greet(LanguageCode lang = LanguageCode::EN) const;
  };

}  // namespace ginger
