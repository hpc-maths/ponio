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
    template <typename solver_t>
    concept has_Snes_method = std::is_member_function_pointer_v<decltype( &solver_t::Snes )>;

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
        solve( operator_t const& op, state_t& u, rhs_t& rhs, std::size_t& n_eval )
        {
            auto solver = samurai::petsc::make_solver( op );
            solver.solve( u, rhs );

            if constexpr ( operator_t::cfg_t::scheme_type == samurai::SchemeType::NonLinear && has_Snes_method<decltype( solver )> )
            {
                int i_n_eval = 0;
                SNESGetNumberFunctionEvals( solver.Snes(), &i_n_eval );

                n_eval = static_cast<std::size_t>( i_n_eval );
            }
            else
            {
                // linear solver, no evaluation of function
                n_eval = 0;
            }

            //::samurai::petsc::solve( op, u, rhs );
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{
    /**
     * @brief For PIROCK method, compute the Shampine's trick
     *
     * @tparam mesh_t  type of underlying mesh in samurai::Field
     * @tparam value_t value stored in samurai::Field
     * @tparam size    the "size" of the samurai::Field (warning: this is not a size in container meaning)
     * @tparam SOA     SOA or AOS
     */
    template <class mesh_t, class value_t, std::size_t size, bool SOA>
    struct shampine_trick<::samurai::Field<mesh_t, value_t, size, SOA>>
    {
        /**
         * @brief solves \f$(I - \alpha R)^{\ell}X = b\f$
         *
         * @tparam l
         * @tparam operator_t
         * @tparam state_t
         * @param alpha           in Shampine's trick \f$\alpha = \gamma \Delta t\f$
         * @param op_reac         \f$R\f$ operator
         * @param initial_guess   initial guess for \f$X\f$ unknown
         * @param rhs             right hand side term, \f$b\f$
         * @param u_tmp           temporary variable
         * @param shampine_result result of unknown \f$X\f$
         */
        template <std::size_t l, typename operator_t, typename state_t>
        void
        operator()( value_t alpha, operator_t&& op_reac, state_t& initial_guess, state_t& rhs, state_t& u_tmp, state_t& shampine_result )
        {
            auto id = ::samurai::make_identity<state_t>();
            // matrix assembly
            auto J_R_op = id - alpha * op_reac;

            auto assembly = samurai::petsc::make_assembly( J_R_op );
            assembly.set_unknown( initial_guess );
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
            Vec rhs_petsc   = samurai::petsc::create_petsc_vector_from( rhs );
            Vec u_tmp_petsc = samurai::petsc::create_petsc_vector_from( result_l1 );
            KSPSolve( ksp, rhs_petsc, u_tmp_petsc );
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
                KSPSolve( ksp, u_tmp_petsc, result_petsc );
            }
        }
    };

} // namespace ponio::shampine_trick
