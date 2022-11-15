#pragma once

#include <fstream>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <filesystem>
#include <type_traits>
#include <iostream>
#include <concepts>
#include <ranges>
#include <iomanip>
#include <limits>


#include "detail.hpp"


namespace observer {

  /** @class capsule
   *  helper to display a value
   *  @tparam state_t   type to display
   *  @details This class only store a reference of data.
   */
  template < typename state_t >
  struct capsule
  {
    state_t const& data;

    capsule( state_t const& dat );
  };

  /**
   * constructor of \ref capsule
   */
  template < typename state_t >
  inline capsule<state_t>::capsule( state_t const& dat ):
    data(dat)
  {}

  /**
   * factory of \ref capsule fron a reference to a data
   */
  template < typename state_t >
  auto
  make_capsule( state_t const& dat ) {
    return capsule< state_t >(dat);   
  }

  #ifndef IN_DOXYGEN
  template < typename state_t >
  std::ostream&
  operator << ( std::ostream & os , capsule<state_t> const& data )
  {
    os << std::setprecision(std::numeric_limits<state_t>::digits10 + 1);
    os << data.data;
    return os;
  }
  #endif

  /**
   * display a \ref capsule which contains an iterable data or not
   */
  template < typename state_t > requires std::ranges::range<state_t>
  std::ostream&
  operator << ( std::ostream & os , capsule<state_t> const& data )
  {
    using value_t = decltype(*std::ranges::cbegin(data.data));

    os << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    std::copy(
      std::ranges::cbegin(data.data), std::ranges::cend(data.data) ,
      std::ostream_iterator<value_t>(os, " ")
    );

    return os;
  }

  /** @class base_observer
   *  base class of observer for CRTP
   *  @tparam Derived_t concret class which store the output stream
   *  @details This class provides a call operator to display in the output stream member of derived class
   */
  template < typename Derived_t >
  struct base_observer
  {
    using derived_t = Derived_t;

    template < typename state_t , typename value_t >
    void
    operator () ( value_t tn , state_t const& un , value_t dt );
  };

  /**
   * call operator to display current state of simulation: `tn un dt` separate only with spaces
   */
  template < typename Derived_t >
  template < typename state_t , typename value_t >
  void
  base_observer<Derived_t>::operator () ( value_t tn , state_t const& un , value_t dt )
  {
    static_cast<derived_t*>(this)->out << tn << " " << make_capsule(un) << " " << dt << "\n";
  }

  /** @class stream_observer
   *  observer that put data into a `std::ostream` done by user
   */
  struct stream_observer : public base_observer<stream_observer>
  {
    std::ostream out;

    stream_observer ( std::ostream & os );
  };

  /**
   * constructor of \ref stream_observer
   * @param os output stream where put the output
   */
  stream_observer::stream_observer ( std::ostream & os ):
    out( os.rdbuf() )
  {}

  /** @class cout_observer
   *  observer that put data into `std::cout`
   */
  struct cout_observer : public stream_observer
  {
    cout_observer ();
  };

  /**
   * constructor of \ref cout_observer
   */
  cout_observer::cout_observer ():
    stream_observer(std::cout)
  {}

  /** @class file_observer
   *  observer that put data into a file
   */
  struct file_observer : public base_observer<file_observer>
  {
    std::ofstream out;

    file_observer ( std::filesystem::path path );

    file_observer ( file_observer const& ) = delete;

    private:

    static std::filesystem::path
    create_directory_if_needed ( std::filesystem::path const& path );
  };

  /**
   * constructor of \ref file_observer
   * @param path path to the output file
   * @note this class creates folder if needed
   */
  file_observer::file_observer ( std::filesystem::path path )
  : out( create_directory_if_needed(path) )
  {}


  /**
   * create parent directory of path if needed and return the path
   * @param path path of output file
   */
  std::filesystem::path
  file_observer::create_directory_if_needed ( std::filesystem::path const& path )
  {
    auto parent = path.parent_path();
    if ( parent.empty() ) {
      std::filesystem::create_directories(parent);
    }
    return path;
  }

  /**
   * litteral to convert a string into \ref file_observer
   */
  file_observer
  operator "" _fobs ( const char * str , std::size_t len )
  {
    return file_observer(std::string_view(str,len));
  }


  /** @class null_observer
   *  observer that do nothing
   */
  struct null_observer
  {
    template < typename state_t , typename value_t >
    void
    operator () ( value_t , state_t const& , value_t );
  };

  /**
   * call operator that do nothing
   */
  template < typename state_t , typename value_t >
  void
  null_observer::operator() ( value_t , state_t const& , value_t )
  {}

  /**
   * @class vector_observer
   * 
   * @tparam state_t type of variable \f$u^n\f$
   * @tparam value_t type of \f$t^n\f$ and \f$\Delta t\f$
   * 
   * This observer saves all iterations inside a vector
   */
  template <typename state_t, typename value_t=double>
  struct vector_observer
  {
    std::vector<std::tuple<value_t,state_t,value_t>> solutions;

    void
    operator () ( value_t tn, state_t const& un, value_t dt );
  };

  /**
   * call operator to save all iteration: `(tn, un, dt)` as a tuple inside a `std::vector`
   */
  template <typename state_t, typename value_t>
  void
  vector_observer<state_t,value_t>::operator () ( value_t tn, state_t const& un, value_t dt )
  {
    solutions.push_back(std::make_tuple(tn,un,dt));
  }

} // namespace observer
