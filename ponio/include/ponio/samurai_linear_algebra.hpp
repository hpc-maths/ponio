// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include "linear_algebra.hpp"

#if __has_include( <samurai/petsc.hpp>)
#include <samurai/petsc.hpp>
#else
#error "Samurai should be included"
#endif

namespace ponio::linear_algebra
{

    template <class mesh_t, class value_t, std::size_t size, bool SOA>
    struct operator_algebra<samurai::Field<mesh_t, value_t, size, SOA>>
    {
        template <typename state_t>
        static auto
        identity( state_t const& )
        {
            return ::samurai::make_identity<state_t>();
        }

        template <typename operator_t, typename state_t, typename rhs_t>
        static void
        solve( operator_t const& op, state_t& u, rhs_t const& rhs )
        {
            ::samurai::petsc::solve( op, u, rhs );
        }
    };

} // namespace ponio::linear_algebra
