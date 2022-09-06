#pragma once

#include <type_traits>

namespace ode {

  /** @class Stage
   *  @brief indicator for overloading each stages of multistages method
   */
  template<std::size_t s>
  using Stage = std::integral_constant<std::size_t, s>;

}
