// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <type_traits>

namespace ponio
{

    /** @class Stage
     *  @brief indicator for overloading each stages of multistages method
     */
    template <std::size_t s>
    using Stage = std::integral_constant<std::size_t, s>;

}
