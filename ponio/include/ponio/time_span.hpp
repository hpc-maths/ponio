// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace ponio
{

    template <typename value_t>
    struct time_span : public std::vector<value_t>
    {
        using std::vector<value_t>::vector;
        using std::vector<value_t>::operator=;

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
