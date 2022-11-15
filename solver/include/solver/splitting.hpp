// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <tuple>
#include <type_traits>
#include <concepts>

namespace ode {
namespace splitting {

  namespace detail {
    template <std::size_t I, typename Problem_t, typename Method_t, typename state_t, typename value_t>
    state_t
    _split_solve( Problem_t & pb, Method_t & meth, state_t & ui, value_t ti, value_t dt )
    {
      value_t current_dt = dt;
      value_t current_time = ti;
      while ( current_time < ti+dt ) {
        std::tie(current_time,ui,current_dt) = std::get<I>(meth)(std::get<I>(pb.system),current_time,ui,current_dt);
        if ( current_time+current_dt > ti+dt ) {
          current_dt = (ti+dt) - current_time;
        }
      }
      return ui;
    }
  }

  /** @class lie
   *  Lie splitting method
   *  @tparam Methods_t list of methods to solve each sub-problem
   */
  template < typename... Methods_t >
  struct lie
  {
    static constexpr std::size_t order = 1;
    static constexpr bool is_splitting_method = true;
    std::tuple<Methods_t...> methods;

    lie ( std::tuple<Methods_t...> const& t );

    /*
    template < std::size_t I=0 , typename Problem_t , typename state_t , typename value_t >
    void _call_step ( Problem_t & f , value_t tn , state_t & ui , value_t dt );
    */

    template < std::size_t I=0 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I == sizeof...(Methods_t)) ,void>::type
    _call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt );
    template < std::size_t I=0 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I < sizeof...(Methods_t)) ,void>::type
    _call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt );

    template < typename Problem_t , typename state_t , typename value_t >
    auto
    operator () ( Problem_t & f , value_t tn , state_t const& un , value_t dt );
  };

  /**
   * constructor of \ref lie from a tuple
   */
  template < typename... Methods_t >
  lie<Methods_t...>::lie ( std::tuple<Methods_t...> const& t ):
    methods(t)
  {}

  /*
   * call current step `I`
   * @tparam I step to call
   * @param f       \ref problem to solve
   * @param tn      current time \f$t^n\f$
   * @param[in,out] ui solution of \f$\dot{u} = \sum_{k=0}^i f_i(t,u)\f$, with initial condition \f$u(t^n)=\texttt{ui}\f$
   * @param dt      time step \f$\Delta t\f$
   * @details compute a solution at time \f$t^n+\Delta t\f$ even with adaptive time step methods
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  void
  lie<Methods_t...>::_call_step ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    value_t current_dt = dt;
    value_t current_time = tn;
    while ( current_time < tn+dt ) {
      std::tie(current_time,ui,current_dt) = std::get<I>(methods)(std::get<I>(f.system),current_time,ui,current_dt);
      if ( current_time+current_dt > tn+dt ) {
        current_dt = (tn+dt) - current_time;
      }
    }
  }
   */

  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I == sizeof...(Methods_t)) ,void>::type
  lie<Methods_t...>::_call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {}

  /**
   * incremental call of each method of each subproblem
   * @tparam I solving step
   * @param f          \ref problem to solve
   * @param tn         current time \f$t^n\f$
   * @param[in,out] ui \f$\texttt{ui}=\phi_{\Delta t}^{[f_1]}\circ\cdots\circ\phi_{\Delta t}^{[f_{i-1}]}(t^n,u^n)\f$
   * @param dt         time step \f$\Delta t\f$
   * @details The parameter @p ui is update to \f$\phi_{\Delta t}^{[f_i]}(t^n,\texttt{ui})\f$
   */
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I < sizeof...(Methods_t)) ,void>::type
  lie<Methods_t...>::_call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    ui = detail::_split_solve<I>(f,methods,ui,tn,dt);
    //_call_step<I>(f,tn,ui,dt);
    _call_inc<I+1>(f,tn,ui,dt);
  }

  /**
   * call operator to initiate Lie splitting recursion
   * @param f  \ref problem to solve
   * @param tn current time \f$t^n\f$
   * @param un current solution \f$u^n \approx u(t^n)\f$ 
   * @param dt time step \f$\Delta t\f$
   */
  template < typename... Methods_t >
  template < typename Problem_t , typename state_t , typename value_t >
  auto
  lie<Methods_t...>::operator () ( Problem_t & f , value_t tn , state_t const& un , value_t dt )
  {
    state_t ui = un;
    _call_inc(f,tn,ui,dt);
    return std::make_tuple(
            tn+dt,
            ui,
            dt
          );
  }

  /**
   * a helper factory for \ref lie functor from a tuple of methods
   * @param t tuple of \ref method
   * @return a \ref lie object build from the tuple of methods
   */
  template < typename... Methods_t >
  auto
  make_lie_from_tuple ( std::tuple<Methods_t...> const& t )
  {
    return lie<Methods_t...>(t);
  }

  /** @class lie_tuple
   *  a helper to deduce method for ::ode::make_method(splitting::lie_tuple<Algorithms_t...> const &, state_t const &)
   *  @tparam Algorithms_t variadic template of algorithms to solve each subproblem
   *  @details This is a dummy class to select correct \ref method to solve the problem
   */
  template < typename... Algorithms_t >
  struct lie_tuple
  {
    static constexpr std::size_t order = 1;
    static constexpr bool is_splitting_method = true;
    std::tuple<Algorithms_t...> algos;

    lie_tuple ( Algorithms_t&&... a );
  };

  /**
   * constructor of \ref lie_tuple from a variadic number of algorithms
   */
  template < typename... Algorithms_t >
  inline lie_tuple<Algorithms_t...>::lie_tuple ( Algorithms_t&&... a ) :
    algos(std::forward<Algorithms_t>(a)...)
  {}

  /**
   * a helper factory for \ref lie_tuple from a tuple of algorithms
   * @param a tuple of \ref method
   * @return a \ref lie_tuple object build from the tuple of methods
   */
  template < typename... Algorithms_t >
  auto
  make_lie_tuple ( Algorithms_t&&... a )
  {
    return lie_tuple<Algorithms_t...>(std::forward<Algorithms_t>(a)...);
  }


  /** @class strang
   *  Strang splitting method
   *  @tparam Methods_t list of methods to solve each sub-problem
   */
  template < typename... Methods_t >
  struct strang : lie<Methods_t...>
  {
    using lie<Methods_t...>::lie;
    using lie<Methods_t...>::methods;
    static constexpr std::size_t order = 2;
    static constexpr bool is_splitting_method = true;

    template < std::size_t I=0 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I == sizeof...(Methods_t)-1) ,void>::type
    _call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt );

    template < std::size_t I=0 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I < sizeof...(Methods_t)-1) ,void>::type
    _call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt );

    template < std::size_t I=sizeof...(Methods_t)-1 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I == 0) ,void>::type
    _call_dec ( Problem_t & f , value_t tn , state_t & ui , value_t dt );

    template < std::size_t I=sizeof...(Methods_t)-1 , typename Problem_t , typename state_t , typename value_t >
    inline typename std::enable_if< (I > 0) ,void>::type
    _call_dec ( Problem_t & f , value_t tn , state_t & ui , value_t dt );

    template < typename Problem_t , typename state_t , typename value_t >
    auto
    operator () ( Problem_t & f , value_t tn , state_t const& un , value_t dt );

  };

  // end of incremental recursion, start of decremental recursion
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I == sizeof...(Methods_t)-1) ,void>::type
  strang<Methods_t...>::_call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    ui = detail::_split_solve<I>(f,methods,ui,tn,dt);
    _call_dec<I-1>(f,tn,ui,dt);
  }

  /**
   * incremental call of each method of each subproblem
   * @tparam I solving step
   * @param f          \ref problem to solve
   * @param tn         current time \f$t^n\f$
   * @param[in,out] ui \f$\texttt{ui}=\phi_{^{\Delta t}/_2}^{[f_1]}\circ\cdots\circ\phi_{^{\Delta t}/_2}^{[f_{i-1}]}(t^n,u^n)\f$
   * @param dt         time step \f$\Delta t\f$
   * @details The parameter @p ui is update to \f$\phi_{^{\Delta t}/_2}^{[f_i]}(t^n,\texttt{ui})\f$
   */
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I < sizeof...(Methods_t)-1) ,void>::type
  strang<Methods_t...>::_call_inc ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    ui = detail::_split_solve<I>(f,methods,ui,tn,0.5*dt);
    _call_inc<I+1>(f,tn,ui,dt);
  }

  // end of decremental recursion, end of recursion
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I == 0) ,void>::type
  strang<Methods_t...>::_call_dec ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    ui = detail::_split_solve<I>(f,methods,ui,tn,0.5*dt);
  }

  /**
   * decremental call of each method of each subproblem
   * @tparam I solving step
   * @param f          \ref problem to solve
   * @param tn         current time \f$t^n\f$
   * @param[in,out] ui \f$\texttt{ui}=\phi_{^{\Delta t}/_2}^{[f_1]}\circ\cdots\circ\phi_{\Delta t}^{[f_{n}]}\circ\cdots\circ\phi_{^{\Delta t}/_2}^{[f_{i+1}]}(t^n,u^n)\f$
   * @param dt         time step \f$\Delta t\f$
   * @details The parameter @p ui is update to \f$\phi_{^{\Delta t}/_2}^{[f_i]}(t^n,\texttt{ui})\f$
   */
  template < typename... Methods_t >
  template < std::size_t I , typename Problem_t , typename state_t , typename value_t >
  inline typename std::enable_if< (I > 0) ,void>::type
  strang<Methods_t...>::_call_dec ( Problem_t & f , value_t tn , state_t & ui , value_t dt )
  {
    ui = detail::_split_solve<I>(f,methods,ui,tn,0.5*dt);
    _call_dec<I-1>(f,tn,ui,dt);
  }

  /**
   * call operator to initiate Strang splitting recursion
   * @param f  \ref problem to solve
   * @param tn current time \f$t^n\f$
   * @param un current solution \f$u^n \approx u(t^n)\f$ 
   * @param dt time step \f$\Delta t\f$
   */
  template < typename... Methods_t >
  template < typename Problem_t , typename state_t , typename value_t >
  auto
  strang<Methods_t...>::operator () ( Problem_t & f , value_t tn , state_t const& un , value_t dt )
  {
    state_t ui = un;
    _call_inc(f,tn,ui,dt);
    return std::make_tuple(
            tn+dt,
            ui,
            dt
          );
  }

  /**
   * a helper factory for \ref strang functor from a tuple of methods
   * @param t tuple of \ref method
   * @return a \ref strang object build from the tuple of methods
   */
  template < typename... Methods_t >
  auto
  make_strang_from_tuple ( std::tuple<Methods_t...> const& t )
  {
    return strang<Methods_t...>(t);
  }

  /** @class strang_tuple
   *  a helper to deduce method for ::ode::make_method(splitting::strang_tuple<Algorithms_t...> const &, state_t const &)
   *  @tparam Algorithms_t variadic template of algorithms to solve each subproblem
   *  @details This is a dummy class to select correct \ref method to solve the problem
   */
  template < typename... Algorithms_t >
  struct strang_tuple
  {
    static constexpr std::size_t order = 2;
    static constexpr bool is_splitting_method = true;
    std::tuple<Algorithms_t...> algos;

    strang_tuple ( Algorithms_t&&... a );
  };

  /**
   * constructor of \ref strang_tuple from a variadic number of algorithms
   */
  template < typename... Algorithms_t >
  inline strang_tuple<Algorithms_t...>::strang_tuple ( Algorithms_t&&... a ) :
    algos(std::forward<Algorithms_t>(a)...)
  {}

  /**
   * a helper factory for \ref strang_tuple from a tuple of algorithms
   * @param a tuple of \ref method
   * @return a \ref strang_tuple object build from the tuple of methods
   */
  template < typename... Algorithms_t >
  auto
  make_strang_tuple ( Algorithms_t&&... a )
  {
    return strang_tuple<Algorithms_t...>(std::forward<Algorithms_t>(a)...);
  }

  template <typename T>
  concept is_splitting_method = requires (T t){ T::is_splitting_method == true; };


} // namespace splitting
} // namespace ode
