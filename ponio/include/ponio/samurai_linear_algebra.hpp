// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cassert>

#include "linear_algebra.hpp"

#if __has_include( <samurai/petsc.hpp>)
#include <xtensor/xfixed.hpp>

#include <samurai/petsc.hpp>
#else
#error "Samurai should be included"
#endif

namespace ponio::linear_algebra
{

    template <class mesh_t, class value_t, std::size_t size, bool SOA>
    struct operator_algebra<::samurai::Field<mesh_t, value_t, size, SOA>>
    {
        template <typename state_t>
        static auto
        identity( state_t const& )
        {
            return ::samurai::make_identity<state_t>();
        }

        template <typename operator_t, typename state_t, typename rhs_t>
        static void
        solve( operator_t const& op, state_t& u, rhs_t& rhs )
        {
            ::samurai::petsc::solve( op, u, rhs );
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{
    template <class mesh_t, class value_t, std::size_t size, bool SOA>
    struct shampine_trick<::samurai::Field<mesh_t, value_t, size, SOA>>
    {
        template <std::size_t l, typename operator_t, typename state_t>
        void
        operator()( value_t gamma_dt, operator_t&& op_reac, state_t& u_sm2pl, state_t& f_D_u, state_t& u_tmp, state_t& shampine_result )
        {
            auto id = ::samurai::make_identity<state_t>();
            // matrix assembly
            auto J_R_op = id - gamma_dt * op_reac;

            auto assembly = samurai::petsc::make_assembly( J_R_op );
            assembly.set_unknown( u_sm2pl );
            Mat J_R;
            assembly.create_matrix( J_R );
            assembly.assemble_matrix( J_R );

            // linear solver
            KSP ksp;
            KSPCreate( PETSC_COMM_SELF, &ksp );
            KSPSetFromOptions( ksp );
            KSPSetOperators( ksp, J_R, J_R );
            PetscInt err = KSPSetUp( ksp );
            if ( err != 0 )
            {
                std::cerr << "The setup of the solver failed!" << std::endl;
                exit( EXIT_FAILURE );
            }

            auto& result_l1 = ( l == 1 ) ? shampine_result : u_tmp;

            // Solve the system
            Vec f_D_U_petsc = samurai::petsc::create_petsc_vector_from( f_D_u );
            Vec u_tmp_petsc = samurai::petsc::create_petsc_vector_from( result_l1 );
            KSPSolve( ksp, u_tmp_petsc, f_D_U_petsc );
            KSPConvergedReason reason_code;
            KSPGetConvergedReason( ksp, &reason_code );
            if ( reason_code < 0 )
            {
                using namespace std::string_literals;
                char const* reason_text;
                KSPGetConvergedReasonString( ksp, &reason_text );
                std::cerr << "Divergence of the solver ("s + reason_text + ")" << std::endl;
                exit( EXIT_FAILURE );
            }

            if constexpr ( l == 2 )
            {
                Vec result_petsc = samurai::petsc::create_petsc_vector_from( shampine_result );
                KSPSolve( ksp, result_petsc, u_tmp_petsc );
            }
        }
    };

} // namespace ponio::shampine_trick
