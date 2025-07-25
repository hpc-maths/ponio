// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cassert>

#include "detail.hpp"
#include "linear_algebra.hpp"

#if __has_include( <samurai/petsc.hpp>)
#include <xtensor/xfixed.hpp>

#include <samurai/petsc.hpp>
#include <samurai/utils.hpp>
#else
#error "Samurai should be included"
#endif

namespace samurai
{
    template <class mesh_t_, class value_t, std::size_t size_, bool SOA>
    class Field;

    template <class mesh_t, class value_t, std::size_t n_comp, bool SOA>
    class VectorField;

    template <class mesh_t, class value_t>
    class ScalarField;
}

namespace ponio_samurai
{
    template <typename field_t>
    concept is_field = std::same_as<field_t,
        ::samurai::Field<typename field_t::mesh_t, typename field_t::value_type, field_t::size, field_t::is_soa>>;

    template <typename field_t>
    concept is_scalar_field = std::same_as<field_t, ::samurai::ScalarField<typename field_t::mesh_t, typename field_t::value_type>>;

    template <typename field_t>
    concept is_vector_field_soa = std::same_as<field_t,
        ::samurai::VectorField<typename field_t::mesh_t, typename field_t::value_type, field_t::n_comp, true>>;
    template <typename field_t>
    concept is_vector_field_aos = std::same_as<field_t,
        ::samurai::VectorField<typename field_t::mesh_t, typename field_t::value_type, field_t::n_comp, false>>;

    template <typename field_t>
    concept is_vector_field = is_vector_field_soa<field_t> || is_vector_field_aos<field_t>;

    template <typename field_t>
    concept is_samurai_field = is_field<field_t> || is_scalar_field<field_t> || is_vector_field<field_t>;
}

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

                VecDestroy( &result_petsc );
            }

            VecDestroy( &u_tmp_petsc );
            VecDestroy( &rhs_petsc );
            KSPDestroy( &ksp );
            MatDestroy( &J_R );
        }
    };

} // namespace ponio::shampine_trick

namespace ponio
{
    /**
     * @brief return a error norm given by: \f$$\sum_i \left(\frac{x_i}{a_{tol} + r_{tol}\max(|y_i|, |z_i|)}\right)^2\Delta x_i\f$$
     *
     * @tparam value_t type of tolerances
     * @tparam state_t type of vectors \f$x\f$, \f$y\f$ and \f$z\f$
     * @param x mainly the estimate error
     * @param y mainly the solution at time \f$t^n\f$
     * @param z mainly the estimation of solution at time \f$t^{n+1}\f$
     * @param a_tol absolute tolerance
     * @param r_tol relative tolerance
     */
    template <typename value_t, typename state_t>
        requires ::ponio_samurai::is_samurai_field<state_t>
    value_t
    norm_error( state_t const& x, state_t const& y, state_t const& z, value_t a_tol, value_t r_tol )
    {
        value_t err = 0.;

        samurai::for_each_cell( x.mesh(),
            [&]( auto& cell )
            {
                err += xt::sum( xt::pow( x[cell] / ( a_tol + r_tol * xt::maximum( xt::abs( y[cell] ), xt::abs( z[cell] ) ) ), 2 ) )[0]
                     * cell.length;
            } );

        return err;
    }

} // namespace ponio
