// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cstddef>

#include <doctest/doctest.h>

#include <ponio/runge_kutta/pirock.hpp>
#include <ponio/runge_kutta/rock_coeff.hpp>

TEST_CASE( "pirock::beta_0_alpha_matches_pirock_f" )
{
    // alpha = 1 / (2 P'_{mdeg+1}(0)) for ell = 1, table recalph of pirock.f (subroutine RECTP),
    // one value per tabulated degree ms[i]; P_{mdeg+1} needs the continuation coefficients recf2
    // clang-format off
    static constexpr std::array<double, 46> recalph = {
        1.0292175779085888, 1.1711426157873857, 1.2891787318896388, 1.3700397217280607,
        1.4302791631568559, 1.4782093408024377, 1.516334027924159,  1.5475988097126931,
        1.5737160337226037, 1.5959315531956058, 1.6150939953884995, 1.6318133984601482,
        1.6465196972300669, 1.659614432610494,  1.67129422100608,   1.6818474251390105,
        1.6914121284472068, 1.7015535221197022, 1.709840773771184,  1.7232927879815139,
        1.7367827079155858, 1.7481922291308314, 1.7580559977582033, 1.7666397849288047,
        1.7741946820706174, 1.7839762257308416, 1.7922570412176893, 1.7993647293953894,
        1.8074256520810972, 1.8142056832103326, 1.819986883764463,  1.826122183655859,
        1.8313031171465501, 1.8357402691831322, 1.8402891266399672, 1.8441678181128414,
        1.8480312856160626, 1.8517598757543596, 1.8552799194744996, 1.8585480291268521,
        1.8615467734561677, 1.864275803993039,  1.866744876869206,  1.869115278174705,
        1.8713339241767621, 1.8732755097069358
    };
    // clang-format on

    using rock_coeff = ponio::runge_kutta::rock::rock2_coeff<double>;
    ponio::runge_kutta::pirock::beta_0<double> const alpha_beta_computer;

    for ( std::size_t i = 0; i < recalph.size(); ++i )
    {
        std::size_t const mdeg = rock_coeff::ms[i];
        std::size_t const s    = mdeg + 2;

        INFO( "mdeg = ", mdeg );
        CHECK( alpha_beta_computer.alpha( s, 1 ) == doctest::Approx( recalph[i] ).epsilon( 1e-12 ) );
    }
}
