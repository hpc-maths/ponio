#pragma once

#include <concepts>

#include "stage.hpp"
#include "detail.hpp"

namespace ode::butcher {

template <typename Tableau>
concept is_embedded_tableau = requires (Tableau t){ t.b2; };

namespace runge_kutta {


  template <typename Tableau>
  struct explicit_rk_butcher
  {
    Tableau butcher;
    static constexpr std::size_t N_stages = Tableau::N_stages;
    static constexpr bool is_embedded = is_embedded_tableau<Tableau>;

    explicit_rk_butcher()
    : butcher()
    {}

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t , std::size_t I >
    inline state_t
    stage ( Stage<I> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return f( tn + butcher.c[I] , ::detail::tpl_inner_product<I>(butcher.A[I], Ki, un, dt) );
    }

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t >
    inline state_t
    stage ( Stage<N_stages> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return ::detail::tpl_inner_product<N_stages>(butcher.b, Ki, un, dt);
    }

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t , typename Tab=Tableau > requires std::same_as<Tab,Tableau> && is_embedded
    inline state_t
    stage ( Stage<N_stages+1> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return ::detail::tpl_inner_product<N_stages>(butcher.b2, Ki, un, dt);
    }
  };

} // namespace runge_kutta

namespace lawson {

  template <typename Exp_t>
  struct lawson_base {
    Exp_t m_exp;

    lawson_base( Exp_t exp_ )
    : m_exp(exp_)
    {}
  };

  template <typename Tableau, typename Exp_t>
  struct explicit_rk_butcher : public lawson_base<Exp_t>
  {
    using lawson_base<Exp_t>::m_exp;

    Tableau butcher;
    static constexpr std::size_t N_stages = Tableau::N_stages;
    static constexpr bool is_embedded = is_embedded_tableau<Tableau>;

    explicit_rk_butcher( Exp_t exp_ )
    : lawson_base<Exp_t>(exp_) , butcher()
    {}

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t , std::size_t i >
    inline state_t
    stage ( Stage<i> , Problem_t & pb , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return m_exp(-butcher.c[i]*dt*pb.l)*pb.n(
        tn + butcher.c[i] ,
        m_exp(butcher.c[i]*dt*pb.l) * ::detail::tpl_inner_product<i>(butcher.A[i], Ki, un, dt)
      );
    }

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t >
    inline state_t
    stage ( Stage<N_stages> , Problem_t & pb , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return m_exp(dt*pb.l) * ::detail::tpl_inner_product<N_stages>(butcher.b, Ki, un, dt);
    }

    template < typename Problem_t , typename state_t , typename value_t , typename ArrayKi_t , typename Tab=Tableau > requires std::same_as<Tab,Tableau> && is_embedded
    inline state_t
    stage ( Stage<N_stages+1> , Problem_t & pb , value_t tn , state_t const& un , ArrayKi_t const& Ki , value_t dt )
    {
      return m_exp(dt*pb.l) * ::detail::tpl_inner_product<N_stages>(butcher.b2, Ki, un, dt);
    }
  };

} // namespace lawson

} // namespace ode::butcher
