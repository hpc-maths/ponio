#pragma once

#include <utility>
#include <functional>
#include <concepts>
#include <array>

namespace detail {

  /* tpl_inner_product */
  template < typename state_t , typename value_t , typename ArrayA_t , typename ArrayB_t , std::size_t ...Is >
  constexpr state_t
  tpl_inner_product_impl ( ArrayA_t const& a , ArrayB_t const& b , state_t const& init , value_t mul_coeff , std::index_sequence<Is...> )
  {
    return ( init + ... + (mul_coeff*a[Is]*b[Is]) );
  }

  /**
   * @brief inner product between two array from 0 to N
   * 
   * @tparam N        number of elements to compute
   * @tparam state_t  type of computed value
   * @tparam value_t  type of coefficents
   * @tparam ArrayA_t type of first array
   * @tparam ArrayB_t type of second array
   * 
   * @param a         first array
   * @param b         second array
   * @param init      starting value to add other values to
   * @param mul_coeff coefficient to multiply each multiplication of inner product
   * 
   * @details This function compute \f$\texttt{init} + \sum_{i=0}^N \texttt{mul_coeff}a_ib_i\f$ without loop thanks to template.
   */
  template < std::size_t N , typename state_t , typename value_t , typename ArrayA_t , typename ArrayB_t >
  constexpr state_t
  tpl_inner_product ( ArrayA_t const& a , ArrayB_t const& b , state_t const& init , value_t mul_coeff=value_t{1.0} )
  {
    return tpl_inner_product_impl(a,b,init,mul_coeff,std::make_index_sequence<N>());
  }


  /* init_fill_array */
  // first version with a value
  template <typename T, std::size_t ...Is>
  constexpr std::array<std::remove_cvref_t<T>, sizeof...(Is)>
  init_fill_array_impl(T && value, std::index_sequence<Is...>)
  {
    return {{ (static_cast<void>(Is), value)... }};
  }
  /**
   * @brief fill an uninitialize array
   * 
   * @tparam T type of stored value in array
   * @tparam N size of array to fill
   * 
   * @param value value to fill in uninitialize array
   * 
   * @code{,cpp}
   *   int i = 42;
   *   const std::array<int,8> arr = detail::init_fill_array<8>( i ); // all values of `arr` are `42`
   * @endcode
   */
  template <std::size_t N, typename T>
  constexpr std::array<std::remove_cvref_t<T>, N>
  init_fill_array(T && value)
  {
    return init_fill_array_impl(std::forward<T>(value), std::make_index_sequence<N>());
  }

  // second version with a invocable parameter (thanks concepts)
  template <typename Function_t, std::size_t ...Is> requires std::invocable<Function_t,std::size_t>
  constexpr std::array<typename decltype(std::function{std::declval<Function_t>()})::result_type, sizeof...(Is)>
  init_fill_array_impl(Function_t && f, std::index_sequence<Is...>)
  {
    return {{ (static_cast<void>(Is), f(Is))... }};
  }
  /**
   * @brief fill an uninitialize array with a generator
   * 
   * @tparam Function_t type of invocable object
   * @tparam N size of array to fill
   * 
   * @param f invocable object that need to get a `std::size_t` and return the value type of output array
   * 
   * @code{,cpp}
   *   const std::array<int,8> arr = detail::init_fill_array<8>([](int i){ return i*i; }); // get {0,1,4,9,16,25,36,49}
   * @endcode
   */
  template <std::size_t N, typename Function_t> requires std::invocable<Function_t,std::size_t>
  constexpr std::array<typename decltype(std::function{std::declval<Function_t>()})::result_type, N>
  init_fill_array(Function_t && f)
  {
    return init_fill_array_impl(std::forward<Function_t>(f), std::make_index_sequence<N>());
  }

  /**
   * @brief concept to check if a template is iterable
   * 
   * @tparam T type of container
   * 
   * @details just check if `std::begin` and `std::end` can be use on type \p T
   */
  template <typename T>
  concept is_iterable = requires (T t) {
    std::begin(t);
    std::end(t);
  };

  /**
   * @brief concept to check if a template is const iterable
   * 
   * @tparam T type of container
   * 
   * @details just check if `std::cbegin` and `std::cend` can be use on type \p T
   */
  template <typename T>
  concept is_const_iterable = requires (T t) {
    std::cbegin(t);
    std::cend(t);
  };


  template <typename Arithmetic, std::size_t ...Is>
  constexpr Arithmetic
  power_impl(Arithmetic && base, std::index_sequence<Is...>)
  {
    return ( static_cast<Arithmetic>(1.0) * ... * (static_cast<void>(Is), base) );
  }

  template <std::size_t Iexp, typename Arithmetic>
  constexpr Arithmetic
  power(Arithmetic && base)
  {
    return power_impl(std::forward<Arithmetic>(base), std::make_index_sequence<Iexp>());
  }

} // namespace detail
