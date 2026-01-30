// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm> // NOLINT(misc-include-cleaner)
#include <cstddef>
#include <vector> // NOLINT(misc-include-cleaner)

namespace ponio
{

    /** @class time_span
     * @brief Class to manage time interval between the start and the end of simulation with optional intermediate values
     *
     * @tparam value_t type of time
     */
    template <typename value_t>
    struct time_span : public std::vector<value_t>
    {
        using std::vector<value_t>::vector;
        using std::vector<value_t>::operator=;

        /**
         * @brief Construct a new time span from a container
         *
         * @tparam Container
         * @param c
         */
        template <typename Container>
        time_span( Container const& c );
    };

    template <typename value_t>
    template <typename Container>
    time_span<value_t>::time_span( Container const& c )
        : std::vector<value_t>( c.size() )
    {
        std::copy( std::begin( c ), std::end( c ), std::vector<value_t>::begin() );
    }

    /**
     * @brief Return evenly spaced numbers over a specified interval
     *
     * @tparam value_t type of time
     * @param start    the stating value of the sequence
     * @param stop     the end value of the sequence
     * @param num      number of samples to generate
     * @param endpoint optional, if true (default) the stop is in the sequence
     */
    template <typename value_t>
    time_span<value_t>
    linspace( value_t start, value_t stop, std::size_t num = 50, bool endpoint = true )
    {
        auto N  = ( endpoint ) ? num - 1 : num;
        auto dt = ( stop - start ) / static_cast<value_t>( N );
        time_span<value_t> t( num );

        for ( std::size_t i = 0; auto& ti : t )
        {
            ti = static_cast<value_t>( i++ ) * dt + start;
        }

        return t;
    }

}
