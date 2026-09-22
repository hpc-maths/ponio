// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cassert>

#include "detail.hpp"
#include "linear_algebra.hpp"

#if __has_include( <samurai/schemes/fv.hpp> )
#include <samurai/field.hpp>
#include <samurai/schemes/fv.hpp>
#include <samurai/utils.hpp>
#else
#error "Samurai should be included"
#endif

namespace ponio_samurai
{
    // Avoid depending on Samurai's internal ScalarField/VectorField template arity.
    template <typename field_t>
    concept is_samurai_field = requires( field_t const& field ) {
                                   typename field_t::mesh_t;
                                   typename field_t::value_type;
                                   field_t::n_comp;
                                   field.mesh();
                               };
}

namespace ponio::detail
{
    /**
     * @brief error estimate of a samurai field
     *
     * The storage of a field spans the ghost cells as well as the cells of the
     * mesh, and the finite volume schemes write the latter only. Walking the
     * storage would therefore weigh stale values in the estimate and average it
     * over more entries than the problem has degrees of freedom, so the cells
     * of the mesh are visited instead.
     *
     * @tparam field_t type of samurai field
     */
    template <typename field_t>
        requires ::ponio_samurai::is_samurai_field<field_t>
    struct error_algebra<field_t>
    {
        template <typename value_t>
        static value_t
        estimate_squared( field_t const& error, field_t const& un, field_t const& unp1, value_t a_tol, value_t r_tol )
        {
            value_t accumulated = static_cast<value_t>( 0 );
            std::size_t n_dof   = 0;

            auto contribution = [&]( value_t e, value_t u_n, value_t u_np1 )
            {
                using namespace std;
                value_t const normalized = abs( e ) / ( a_tol + r_tol * max( abs( u_n ), abs( u_np1 ) ) );

                accumulated += normalized * normalized;
                ++n_dof;
            };

            ::samurai::for_each_cell( error.mesh(),
                [&]( auto& cell )
                {
                    if constexpr ( field_t::n_comp == 1 )
                    {
                        contribution( error[cell], un[cell], unp1[cell] );
                    }
                    else
                    {
                        for ( std::size_t c = 0; c < field_t::n_comp; ++c )
                        {
                            contribution( error[cell]( c ), un[cell]( c ), unp1[cell]( c ) );
                        }
                    }
                } );

            return accumulated / static_cast<value_t>( n_dof );
        }
    };

} // namespace ponio::detail

namespace ponio::linear_algebra
{
    template <typename solver_t>
    concept has_Snes_method = std::is_member_function_pointer_v<decltype( &solver_t::Snes )>;

    template <typename field_t>
        requires ::ponio_samurai::is_samurai_field<field_t>
    struct operator_algebra<field_t>
    {
        template <typename state_t>
        static auto
        identity( state_t const& )
        {
            return ::samurai::make_identity<state_t>();
        }

        template <typename operator_t, typename state_t, typename rhs_t>
        static void
        solve( operator_t& op, state_t& u, rhs_t& rhs, std::size_t& n_eval )
        {
            auto solver = samurai::petsc::make_solver( op );
            auto _rhs   = static_cast<state_t>( rhs );
            solver.solve( u, _rhs );

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

    /**
     * @brief in-place update of a samurai field
     *
     * samurai fields provide addition and assignment but no compound
     * assignment, so both operations are written as an assignment from an
     * expression. Assignment is elementwise, so the field appearing on both
     * sides reads and writes the same cell and needs no temporary.
     *
     * @tparam field_t type of samurai field
     */
    template <typename field_t>
        requires ::ponio_samurai::is_samurai_field<field_t>
    struct state_algebra<field_t>
    {
        template <typename value_t>
        static void
        scale( field_t& y, value_t alpha )
        {
            y = alpha * y;
        }

        template <typename increment_t>
        static void
        add( field_t& y, increment_t const& x )
        {
            y = y + x;
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{
    /**
     * @brief For PIROCK method, compute the Shampine's trick
     *
     * @tparam field_t type of samurai field
     */
    template <typename field_t>
        requires ::ponio_samurai::is_samurai_field<field_t>
    struct shampine_trick<field_t>
    {
        using value_t = typename field_t::value_type;

        /**
         * @brief solves \f$(I - \alpha R)^{\ell}X = b\f$
         *
         * @tparam ell
         * @tparam operator_t
         * @tparam state_t
         * @param alpha           in Shampine's trick \f$\alpha = \gamma \Delta t\f$
         * @param op_reac         \f$R\f$ operator
         * @param initial_guess   initial guess for \f$X\f$ unknown
         * @param rhs             right hand side term, \f$b\f$
         * @param u_tmp           temporary variable
         * @param shampine_result result of unknown \f$X\f$
         */
        template <std::size_t ell, typename operator_t, typename state_t>
        void
        operator()( value_t alpha, operator_t&& op_reac, state_t& initial_guess, state_t& rhs, state_t& u_tmp, state_t& shampine_result )
        {
            bool constexpr is_local_operator = operator_t::cfg_t::stencil_size == 1;

            auto id = ::samurai::make_identity<state_t>();
            // matrix assembly
            auto J_R_op = id - alpha * op_reac;

            auto assembly = samurai::petsc::make_assembly( J_R_op );
            assembly.set_unknown( initial_guess );

            if ( is_local_operator )
            {
                assembly.include_bc( false );
            }

            Mat J_R;
            assembly.create_matrix( J_R );
            assembly.assemble_matrix( J_R );

            // linear solver
            KSP ksp;
            KSPCreate( PETSC_COMM_SELF, &ksp );
            KSPSetFromOptions( ksp );
            KSPSetOperators( ksp, J_R, J_R );
            PetscInt const err = KSPSetUp( ksp );
            if ( err != 0 )
            {
                std::cerr << "The setup of the solver failed!" << std::endl;
                exit( EXIT_FAILURE ); // NOLINT
            }

            auto& result_l1 = ( ell == 1 ) ? shampine_result : u_tmp;

            // Solve the system
            Vec rhs_petsc   = samurai::petsc::create_petsc_vector_from( rhs );
            Vec u_tmp_petsc = samurai::petsc::create_petsc_vector_from( result_l1 );

            assembly.set_0_for_all_ghosts( rhs_petsc );

            KSPSolve( ksp, rhs_petsc, u_tmp_petsc );
            KSPConvergedReason reason_code;
            KSPGetConvergedReason( ksp, &reason_code );
            if ( reason_code < 0 )
            {
                using namespace std::string_literals;
                char const* reason_text;
                KSPGetConvergedReasonString( ksp, &reason_text );
                std::cerr << "Divergence of the solver ("s + reason_text + ")" << std::endl;
                exit( EXIT_FAILURE ); // NOLINT
            }

            if constexpr ( ell == 2 )
            {
                Vec result_petsc = samurai::petsc::create_petsc_vector_from( shampine_result );

                assembly.set_0_for_all_ghosts( u_tmp_petsc );

                KSPSolve( ksp, u_tmp_petsc, result_petsc );
                KSPGetConvergedReason( ksp, &reason_code );
                if ( reason_code < 0 )
                {
                    using namespace std::string_literals;
                    char const* reason_text;
                    KSPGetConvergedReasonString( ksp, &reason_text );
                    std::cerr << "Divergence of the solver ("s + reason_text + ")" << std::endl;
                    exit( EXIT_FAILURE ); // NOLINT
                }

                VecDestroy( &result_petsc );
            }

            VecDestroy( &u_tmp_petsc );
            VecDestroy( &rhs_petsc );
            KSPDestroy( &ksp );
            MatDestroy( &J_R );
        }
    };

} // namespace ponio::shampine_trick
