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
    static constexpr std::size_t order = Tableau::order;
    static constexpr const char* id = Tableau::id;

    explicit_rk_butcher(double tol_=1e-4)
    : butcher(), tol(tol_)
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

    double tol;
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

  /**
   * factory of Lawson method to help clang++ which doesn't support C++20
   * 
   * @tparam Tableau type of Butcher tableau
   * @tparam Exp_t type of exponential function
   * @param exp_ exponential function
   */
  template <typename Tableau, typename Exp_t>
  auto
  make_lawson( Exp_t exp_ )
  {
    return explicit_rk_butcher<Tableau,Exp_t>(exp_);
  }

} // namespace lawson

namespace chebyshev {

  /**
   * @class T
   * Chebyshev polynomial of first kind by recursive method: \f$$\begin{aligned}T_0(x) &= 1\\ T_1(x) &= x\\ T_{n+1}(x) &= 2xT_n(x) - T_{n-1}(x)\end{aligned}}\f$$
   * 
   * @tparam N degree of polynomial
   */
  template <std::size_t N>
  struct T
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x );
  };

  /**
   * value of Chebyshev polynomial \f$T_N(x)\f$
   * 
   * @tparam value_t type of \f$x\f$
   * @param x value where evaluate \f$T_N\f$
   */
  template <std::size_t N>
  template <typename value_t>
  constexpr value_t
  T<N>::value ( value_t x )
  {
    return static_cast<value_t>(2.)*x*T<N-1>::value(x) - T<N-2>::value(x);
  }

  template <>
  struct T<0>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return static_cast<value_t>(1.);
    }
  };

  template <>
  struct T<1>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return x;
    }
  };

  /**
   * @class U
   * Chebyshev polynomial of second kind by recursive method: \f$$\begin{aligned}U_0(x) &= 1\\ U_1(x) &= 2x\\ U_{n+1}(x) &= 2xU_n(x) - U_{n-1}(x)\end{aligned}}\f$$
   * 
   * @tparam N degree of polynomial
   */
  template <std::size_t N>
  struct U
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x );
  };

  /**
   * value of Chebyshev polynomial \f$U_N(x)\f$
   * 
   * @tparam value_t type of \f$x\f$
   * @param x value where evaluate \f$U_N\f$
   */
  template <std::size_t N>
  template <typename value_t>
  constexpr value_t
  U<N>::value ( value_t x )
  {
    return static_cast<value_t>(2.)*x*U<N-1>::value(x) - U<N-2>::value(x);
  }

  template <>
  struct U<0>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return static_cast<value_t>(1.);
    }
  };

  template <>
  struct U<1>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return static_cast<value_t>(2.)*x;
    }
  };

  /**
   * @class dT
   * derivatives of Chebyshev polynomial: \f$\frac{\mathrm{d}T_N}{\mathrm{d}x}(x)=NU_{N-1}(x)\f$
   * 
   * @tparam N degree of polynomial
   */
  template <std::size_t N>
  struct dT
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x );
  };

  /**
   * value of derivatives of Chebyshev polynomial: \f$\frac{\mathrm{d}T_N}{\mathrm{d}x}(x)\f$
   * 
   * @tparam value_t type of \f$x\f$
   * @param x value where evaluate \f$\frac{\mathrm{d}T_N}{\mathrm{d}x}\f$
   */
  template <std::size_t N>
  template <typename value_t>
  constexpr value_t
  dT<N>::value ( value_t x )
  {
    return N*U<N-1>::value(x);
  }

  template <>
  struct dT<0>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return static_cast<value_t>(0.);
    }
  };

  /**
   * @class ddT
   * second derivative of Chebyshev polynomial: \f$\frac{\mathrm{d}^2T_N}{\mathrm{d}x^2}(x)=N\frac{NT_N(x) - xU_{N-1}(x)}{x^2 -1}\f$
   * 
   * @tparam N degree of polynomial
   */
  template <std::size_t N>
  struct ddT
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x );
  };

  /**
   * value of derivative of Chebyshev polynomial: \f$\frac{\mathrm{d}T_N}{\mathrm{d}x}(x)\f$
   * 
   * @tparam value_t type of \f$x\f$
   * @param x value where evaluate \f$\frac{\mathrm{d}T_N}{\mathrm{d}x}\f$
   */
  template <std::size_t N>
  template <typename value_t>
  constexpr value_t
  ddT<N>::value ( value_t x )
  {
    return N*(N*T<N>::value(x) - x*U<N-1>::value(x))/(x*x - static_cast<value_t>(1.));
  }

  template <>
  struct ddT<0>
  {
    template <typename value_t>
    static constexpr value_t
    value ( value_t x )
    {
      return static_cast<value_t>(0.);
    }
  };



  /** @class explicit_rkc2
   *  @brief define RKC2 with user defined number of stages
   * 
   *  @tparam _Nstages number of stages
   *  @tparam value_t type of coefficients
   */
  template <std::size_t _N_stages, typename value_t=double>
  struct explicit_rkc2
  {
    static constexpr bool is_embedded = false;
    static constexpr std::size_t N_stages = _N_stages;
    value_t w0;
    value_t w1;

    template <std::size_t J, typename T=value_t>
    struct b
    {
      static constexpr value_t
      value( value_t x )
      {
        return ddT<J>::value(x) / ( dT<J>::value(x)*dT<J>::value(x) );
      }
    };
    template <typename T>
    struct b<0,T>
    {
      static constexpr value_t
      value( value_t x )
      {
        return b<2,value_t>::value(x);
      }
    };
    template <typename T>
    struct b<1,T>
    {
      static constexpr value_t
      value( value_t x )
      {
        return b<2,value_t>::value(x);
      }
    };


    explicit_rkc2 ( value_t eps=static_cast<value_t>(2./13.) )
    : w0(1. + eps/(N_stages*N_stages))
    {
      w1 = dT<N_stages>::value(w0)/ddT<N_stages>::value(w0);
    }

    template < typename Problem_t , typename state_t , typename ArrayKi_t , std::size_t j >
    inline state_t
    stage ( Stage<j> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Yi , value_t dt )
    {
      value_t mj  = 2.*b<j>::value(w0)/b<j-1>::value(w0)*w0;
      value_t nj  =   -b<j>::value(w0)/b<j-2>::value(w0);
      value_t mjt = 2.*b<j>::value(w0)/b<j-1>::value(w0)*w1;
      value_t gjt = -(1. - b<j-1>::value(w0)*T<j-1>::value(w0))*mjt;
      value_t cjm1 = dT<N_stages>::value(w0)/ddT<N_stages>::value(w0)*ddT<j-1>::value(w0)/dT<j-1>::value(w0);

      return (1.-mj-nj)*un + mj*Yi[j-1] + nj*Yi[j-2] + mjt*dt*f(tn+cjm1*dt,Yi[j-1]) + gjt*dt*Yi[0];
    }

    template < typename Problem_t , typename state_t , typename ArrayKi_t >
    inline state_t
    stage ( Stage<0> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Yi , value_t dt )
    {
      return f(tn,un); // be careful Yi[0] stores f(tn,un) not un!!!
    }

    template < typename Problem_t , typename state_t , typename ArrayKi_t >
    inline state_t
    stage ( Stage<1> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Yi , value_t dt )
    {
      value_t m1t = b<1>::value(w0)*w1;
      return un + dt*m1t*Yi[0];
    }

    template < typename Problem_t , typename state_t , typename ArrayKi_t >
    inline state_t
    stage ( Stage<2> , Problem_t & f , value_t tn , state_t const& un , ArrayKi_t const& Yi , value_t dt )
    {
      value_t m2 = 2.*w0;
      value_t n2 = -1.;
      value_t m2t = 2.*w1;
      value_t c2 = dT<N_stages>::value(w0)/ddT<N_stages>::value(w0)*ddT<2>::value(w0)/dT<2>::value(w0);
      value_t c1 = c2/dT<2>::value(w0);
      value_t g2t = -(1. - b<1>::value(w0)*T<1>::value(w0))*m2t;

      return (1.-m2-n2)*un + m2*Yi[1] + n2*un + m2t*dt*f(tn+c1*dt,Yi[1]) + g2t*dt*Yi[0];
    }
  };

} // namespace chebyshev

} // namespace ode::butcher
