#pragma once

#include "butcher_tableau.hpp"
#include "generic_butcher_rk.hpp"

namespace ode::butcher {


template <typename value_t=double>
struct butcher_euler : public butcher_tableau<1,value_t>
{
  using base_t = butcher_tableau<1,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_euler()
  : base_t(
    {{
        { 0 }
    }}, // A
    { 1 }, // b
    { 0 }  // c
  )
  {}
};
template <typename value_t=double>
using euler = runge_kutta::explicit_rk_butcher<butcher_euler<value_t>>;


template <typename value_t=double>
struct butcher_explicit_euler_sub4 : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_explicit_euler_sub4()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/4 , 0 , 0 , 0 },{ 1/4 , 1/4 , 0 , 0 },{ 1/4 , 1/4 , 1/4 , 0 }
    }}, // A
    { 1/4 , 1/4 , 1/4 , 1/4 }, // b
    { 0 , 1/4 , 1/2 , 3/4 }  // c
  )
  {}
};
template <typename value_t=double>
using explicit_euler_sub4 = runge_kutta::explicit_rk_butcher<butcher_explicit_euler_sub4<value_t>>;


template <typename value_t=double>
struct butcher_rk54_6m : public butcher_tableau<6,value_t>
{
  using base_t = butcher_tableau<6,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk54_6m()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 },{ 1/5 , 0 , 0 , 0 , 0 , 0 },{ 3/40 , 9/40 , 0 , 0 , 0 , 0 },{ 3/10 , -9/10 , 6/5 , 0 , 0 , 0 },{ 226/729 , -25/27 , 880/729 , 55/729 , 0 , 0 },{ -181/270 , 5/2 , -266/297 , -91/27 , 189/55 , 0 }
    }}, // A
    { 19/216 , 0 , 1000/2079 , -125/216 , 81/88 , 5/56 }, // b
    { 0 , 1/5 , 3/10 , 3/5 , 2/3 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk54_6m = runge_kutta::explicit_rk_butcher<butcher_rk54_6m<value_t>>;


template <typename value_t=double>
struct butcher_rk54_7m : public butcher_tableau<7,value_t>
{
  using base_t = butcher_tableau<7,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk54_7m()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/5 , 0 , 0 , 0 , 0 , 0 , 0 },{ 3/40 , 9/40 , 0 , 0 , 0 , 0 , 0 },{ 44/45 , -56/15 , 32/9 , 0 , 0 , 0 , 0 },{ 19372/6561 , -25360/2187 , 64448/6561 , -212/729 , 0 , 0 , 0 },{ 9017/3168 , -355/33 , 46732/5247 , 49/176 , -5103/18656 , 0 , 0 },{ 35/384 , 0 , 500/1113 , 125/192 , -2187/6784 , 11/84 , 0 }
    }}, // A
    { 35/384 , 0 , 500/1113 , 125/192 , -2187/6784 , 11/84 , 0 }, // b
    { 0 , 1/5 , 3/10 , 4/5 , 8/9 , 1 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk54_7m = runge_kutta::explicit_rk_butcher<butcher_rk54_7m<value_t>>;


template <typename value_t=double>
struct butcher_rk54_7s : public butcher_tableau<7,value_t>
{
  using base_t = butcher_tableau<7,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk54_7s()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 2/9 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/12 , 1/4 , 0 , 0 , 0 , 0 , 0 },{ 55/324 , -25/108 , 50/81 , 0 , 0 , 0 , 0 },{ 83/330 , -13/22 , 61/66 , 9/110 , 0 , 0 , 0 },{ -19/28 , 9/4 , 1/7 , -27/7 , 22/7 , 0 , 0 },{ 19/200 , 0 , 3/5 , -243/400 , 33/40 , 7/80 , 0 }
    }}, // A
    { 19/200 , 0 , 3/5 , -243/400 , 33/40 , 7/80 , 0 }, // b
    { 0 , 2/9 , 1/3 , 5/9 , 2/3 , 1 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk54_7s = runge_kutta::explicit_rk_butcher<butcher_rk54_7s<value_t>>;


template <typename value_t=double>
struct butcher_rk6es : public butcher_tableau<7,value_t>
{
  using base_t = butcher_tableau<7,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk6es()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 0.202276644898141 , 0 , 0 , 0 , 0 , 0 , 0 },{ 0.0758537418368027 , 0.227561225510408 , 0 , 0 , 0 , 0 , 0 },{ 1.35928221728330 , -5.23788570262881 , 4.75360348534551 , 0 , 0 , 0 , 0 },{ -0.321092002258022 , 1.65135312792238 , -0.905286676763721 , 0.0750255510993598 , 0 , 0 , 0 },{ 0.292321839349364 , -0.748269386089830 , 0.592470844966485 , -0.0395545388491436 , 0.0280312406231245 , 0 , 0 },{ -20.6627618949041 , 63.8523209463321 , -74.1517509476888 , 0.864117644373384 , 14.5054816592948 , 16.5925925925926 , 0 }
    }}, // A
    { 0.0142857142857143 , 0 , 0 , 0.270899470899471 , 0.429629629629630 , 0.270899470899471 , 0.0142857142857143 }, // b
    { 0 , 0.202276644898141 , 0.303414967347211 , 0.875000000000000 , 0.500000000000000 , 0.125000000000000 , 1.00000000000000 }  // c
  )
  {}
};
template <typename value_t=double>
using rk6es = runge_kutta::explicit_rk_butcher<butcher_rk6es<value_t>>;


template <typename value_t=double>
struct butcher_rk_118 : public butcher_tableau<11,value_t>
{
  using base_t = butcher_tableau<11,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_118()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/4 , 1/4 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/7 , 3*std::pow(21, 1/2)/98 - 1/14 , 3/7 - 5*std::pow(21, 1/2)/49 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 11/84 - std::pow(21, 1/2)/84 , 0 , 2/7 - 4*std::pow(21, 1/2)/63 , std::pow(21, 1/2)/252 + 1/12 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 5/48 - std::pow(21, 1/2)/48 , 0 , 1/4 - std::pow(21, 1/2)/36 , -7*std::pow(21, 1/2)/180 - 77/120 , 7*std::pow(21, 1/2)/80 + 63/80 , 0 , 0 , 0 , 0 , 0 , 0 },{ std::pow(21, 1/2)/42 + 5/21 , 0 , -92*std::pow(21, 1/2)/315 - 48/35 , 29*std::pow(21, 1/2)/18 + 211/30 , -23*std::pow(21, 1/2)/14 - 36/5 , 13*std::pow(21, 1/2)/35 + 9/5 , 0 , 0 , 0 , 0 , 0 },{ 1/14 , 0 , 0 , 0 , std::pow(21, 1/2)/42 + 1/9 , std::pow(21, 1/2)/21 + 13/63 , 1/9 , 0 , 0 , 0 , 0 },{ 1/32 , 0 , 0 , 0 , 7*std::pow(21, 1/2)/192 + 91/576 , 11/72 , 25*std::pow(21, 1/2)/384 - 385/1152 , 63/128 - 13*std::pow(21, 1/2)/128 , 0 , 0 , 0 },{ 1/14 , 0 , 0 , 0 , 1/9 , std::pow(21, 1/2)/15 - 733/2205 , 515/504 - 37*std::pow(21, 1/2)/168 , 11*std::pow(21, 1/2)/56 - 51/56 , 132/245 - 4*std::pow(21, 1/2)/35 , 0 , 0 },{ 0 , 0 , 0 , 0 , -7*std::pow(21, 1/2)/18 - 7/3 , -28*std::pow(21, 1/2)/45 - 2/5 , 53*std::pow(21, 1/2)/72 - 91/24 , 301/72 - 53*std::pow(21, 1/2)/72 , 28*std::pow(21, 1/2)/45 + 28/45 , 7*std::pow(21, 1/2)/18 + 49/18 , 0 }
    }}, // A
    { 1/20 , 0 , 0 , 0 , 0 , 0 , 0 , 49/180 , 16/45 , 49/180 , 1/20 }, // b
    { 0 , 1/2 , 1/2 , 1/2 - std::pow(21, 1/2)/14 , 1/2 - std::pow(21, 1/2)/14 , 1/2 , std::pow(21, 1/2)/14 + 1/2 , std::pow(21, 1/2)/14 + 1/2 , 1/2 , 1/2 - std::pow(21, 1/2)/14 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_118 = runge_kutta::explicit_rk_butcher<butcher_rk_118<value_t>>;


template <typename value_t=double>
struct butcher_rk_21 : public butcher_tableau<2,value_t>
{
  using base_t = butcher_tableau<2,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_21()
  : base_t(
    {{
        { 0 , 0 },{ 1 , 0 }
    }}, // A
    { 1/2 , 1/2 }, // b
    { 0 , 1/2 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_21 = runge_kutta::explicit_rk_butcher<butcher_rk_21<value_t>>;


template <typename value_t=double>
struct butcher_rk_22_midpoint : public butcher_tableau<2,value_t>
{
  using base_t = butcher_tableau<2,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_22_midpoint()
  : base_t(
    {{
        { 0 , 0 },{ 1/2 , 0 }
    }}, // A
    { 0 , 1 }, // b
    { 0 , 1/2 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_22_midpoint = runge_kutta::explicit_rk_butcher<butcher_rk_22_midpoint<value_t>>;


template <typename value_t=double>
struct butcher_rk_22_ralston : public butcher_tableau<2,value_t>
{
  using base_t = butcher_tableau<2,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_22_ralston()
  : base_t(
    {{
        { 0 , 0 },{ 2/3 , 0 }
    }}, // A
    { 1/4 , 3/4 }, // b
    { 0 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_22_ralston = runge_kutta::explicit_rk_butcher<butcher_rk_22_ralston<value_t>>;


template <typename value_t=double>
struct butcher_rk_32_best : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_32_best()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/2 , 0 , 0 },{ 0 , 1/2 , 0 }
    }}, // A
    { 0 , 0 , 1 }, // b
    { 0 , 1/2 , 1/2 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_32_best = runge_kutta::explicit_rk_butcher<butcher_rk_32_best<value_t>>;


template <typename value_t=double>
struct butcher_rk_33 : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/2 , 0 , 0 },{ -1 , 2 , 0 }
    }}, // A
    { 1/6 , 2/3 , 1/6 }, // b
    { 0 , 1/2 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33 = runge_kutta::explicit_rk_butcher<butcher_rk_33<value_t>>;


template <typename value_t=double>
struct butcher_rk_33_233e : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33_233e()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 2/3 , 0 , 0 },{ 1/3 , 1/3 , 0 }
    }}, // A
    { 1/4 , 0 , 3/4 }, // b
    { 0 , 2/3 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33_233e = runge_kutta::explicit_rk_butcher<butcher_rk_33_233e<value_t>>;


template <typename value_t=double>
struct butcher_rk_33_bogackishampine : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33_bogackishampine()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/2 , 0 , 0 },{ 0 , 3/4 , 0 }
    }}, // A
    { 2/9 , 1/3 , 4/9 }, // b
    { 0 , 1/2 , 3/4 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33_bogackishampine = runge_kutta::explicit_rk_butcher<butcher_rk_33_bogackishampine<value_t>>;


template <typename value_t=double>
struct butcher_rk_33_heun : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33_heun()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/3 , 0 , 0 },{ 0 , 2/3 , 0 }
    }}, // A
    { 1/4 , 0 , 3/4 }, // b
    { 0 , 1/3 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33_heun = runge_kutta::explicit_rk_butcher<butcher_rk_33_heun<value_t>>;


template <typename value_t=double>
struct butcher_rk_33_ralston : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33_ralston()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/2 , 0 , 0 },{ 0 , 3/4 , 0 }
    }}, // A
    { 2/9 , 1/3 , 4/9 }, // b
    { 0 , 1/2 , 3/4 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33_ralston = runge_kutta::explicit_rk_butcher<butcher_rk_33_ralston<value_t>>;


template <typename value_t=double>
struct butcher_rk_33_van_der_houwen : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_33_van_der_houwen()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 8/15 , 0 , 0 },{ 1/4 , 5/12 , 0 }
    }}, // A
    { 1/4 , 0 , 3/4 }, // b
    { 0 , 8/15 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_33_van_der_houwen = runge_kutta::explicit_rk_butcher<butcher_rk_33_van_der_houwen<value_t>>;


template <typename value_t=double>
struct butcher_rk_44 : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_44()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/2 , 0 , 0 , 0 },{ 0 , 1/2 , 0 , 0 },{ 0 , 0 , 1 , 0 }
    }}, // A
    { 1/6 , 1/3 , 1/3 , 1/6 }, // b
    { 0 , 1/2 , 1/2 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_44 = runge_kutta::explicit_rk_butcher<butcher_rk_44<value_t>>;


template <typename value_t=double>
struct butcher_rk_44_235j : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_44_235j()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/4 , 0 , 0 , 0 },{ 0 , 1/2 , 0 , 0 },{ 1 , -2 , 2 , 0 }
    }}, // A
    { 1/6 , 0 , 2/3 , 1/6 }, // b
    { 0 , 1/4 , 1/2 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_44_235j = runge_kutta::explicit_rk_butcher<butcher_rk_44_235j<value_t>>;


template <typename value_t=double>
struct butcher_rk_44_38 : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_44_38()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/3 , 0 , 0 , 0 },{ -1/3 , 1 , 0 , 0 },{ 1 , -1 , 1 , 0 }
    }}, // A
    { 1/8 , 3/8 , 3/8 , 1/8 }, // b
    { 0 , 1/3 , 2/3 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_44_38 = runge_kutta::explicit_rk_butcher<butcher_rk_44_38<value_t>>;


template <typename value_t=double>
struct butcher_rk_44_ralston : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_44_ralston()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 0.400000000000000 , 0 , 0 , 0 },{ 0.296977610000000 , 0.158759640000000 , 0 , 0 },{ 0.218100400000000 , -3.05096516000000 , 3.83286476000000 , 0 }
    }}, // A
    { 0.174760280000000 , -0.551480660000000 , 1.20553560000000 , 0.171184780000000 }, // b
    { 0 , 0.400000000000000 , -0.455737250000000 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_44_ralston = runge_kutta::explicit_rk_butcher<butcher_rk_44_ralston<value_t>>;


template <typename value_t=double>
struct butcher_rk_65 : public butcher_tableau<6,value_t>
{
  using base_t = butcher_tableau<6,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_65()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 },{ 1/5 , 0 , 0 , 0 , 0 , 0 },{ 3/40 , 9/40 , 0 , 0 , 0 , 0 },{ 44/45 , -56/15 , 32/9 , 0 , 0 , 0 },{ 19372/6561 , -25360/2187 , 64448/6561 , -212/729 , 0 , 0 },{ 9017/3168 , -355/33 , 46732/5247 , 49/176 , -5103/18656 , 0 }
    }}, // A
    { 35/384 , 0 , 500/1113 , 125/192 , -2187/6784 , 11/84 }, // b
    { 0 , 1/5 , 3/10 , 4/5 , 8/9 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_65 = runge_kutta::explicit_rk_butcher<butcher_rk_65<value_t>>;


template <typename value_t=double>
struct butcher_rk_65_236a : public butcher_tableau<6,value_t>
{
  using base_t = butcher_tableau<6,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_65_236a()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 },{ 1/4 , 0 , 0 , 0 , 0 , 0 },{ 1/8 , 1/8 , 0 , 0 , 0 , 0 },{ 0 , 0 , 1/2 , 0 , 0 , 0 },{ 3/16 , -3/8 , 3/8 , 9/16 , 0 , 0 },{ -3/7 , 8/7 , 6/7 , -12/7 , 8/7 , 0 }
    }}, // A
    { 7/90 , 0 , 16/45 , 2/15 , 16/45 , 7/90 }, // b
    { 0 , 1/4 , 1/4 , 1/2 , 3/4 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_65_236a = runge_kutta::explicit_rk_butcher<butcher_rk_65_236a<value_t>>;


template <typename value_t=double>
struct butcher_rk_76 : public butcher_tableau<7,value_t>
{
  using base_t = butcher_tableau<7,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_76()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1 , 0 , 0 , 0 , 0 , 0 , 0 },{ 3/8 , 1/8 , 0 , 0 , 0 , 0 , 0 },{ 8/27 , 2/27 , 8/27 , 0 , 0 , 0 , 0 },{ 9*std::pow(21, 1/2)/392 - 3/56 , std::pow(21, 1/2)/49 - 1/7 , 6/7 - 6*std::pow(21, 1/2)/49 , 3*std::pow(21, 1/2)/392 - 9/56 , 0 , 0 , 0 },{ -51*std::pow(21, 1/2)/392 - 33/56 , -std::pow(21, 1/2)/49 - 1/7 , -8*std::pow(21, 1/2)/49 , 363*std::pow(21, 1/2)/1960 + 9/280 , std::pow(21, 1/2)/5 + 6/5 , 0 , 0 },{ 7*std::pow(21, 1/2)/12 + 11/6 , 2/3 , 14*std::pow(21, 1/2)/9 - 10/9 , 7/10 - 21*std::pow(21, 1/2)/20 , -7*std::pow(21, 1/2)/10 - 343/90 , 49/18 - 7*std::pow(21, 1/2)/18 , 0 }
    }}, // A
    { 1/20 , 0 , 16/45 , 0 , 49/180 , 49/180 , 1/20 }, // b
    { 0 , 1 , 1/2 , 2/3 , 1/2 - std::pow(21, 1/2)/14 , std::pow(21, 1/2)/14 + 1/2 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_76 = runge_kutta::explicit_rk_butcher<butcher_rk_76<value_t>>;


template <typename value_t=double>
struct butcher_rk_86 : public butcher_tableau<8,value_t>
{
  using base_t = butcher_tableau<8,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_86()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/9 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/24 , 1/8 , 0 , 0 , 0 , 0 , 0 , 0 },{ 1/6 , -1/2 , 2/3 , 0 , 0 , 0 , 0 , 0 },{ 935/2536 , -2781/2536 , 309/317 , 321/1268 , 0 , 0 , 0 , 0 },{ -12710/951 , 8287/317 , -40/317 , -6335/317 , 8 , 0 , 0 , 0 },{ 5840285/3104064 , -7019/2536 , -52213/86224 , 1278709/517344 , -433/2448 , 33/1088 , 0 , 0 },{ -5101675/1767592 , 112077/25994 , 334875/441898 , -973617/883796 , -1421/1394 , 333/5576 , 36/41 , 0 }
    }}, // A
    { 41/840 , 0 , 9/35 , 9/280 , 34/105 , 9/280 , 9/35 , 41/840 }, // b
    { 0 , 1/9 , 1/6 , 1/3 , 1/2 , 2/3 , 5/6 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_86 = runge_kutta::explicit_rk_butcher<butcher_rk_86<value_t>>;


template <typename value_t=double>
struct butcher_rk_nssp_21 : public butcher_tableau<2,value_t>
{
  using base_t = butcher_tableau<2,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_nssp_21()
  : base_t(
    {{
        { 0 , 0 },{ 3/4 , 0 }
    }}, // A
    { 0 , 1 }, // b
    { 0 , 3/4 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_nssp_21 = runge_kutta::explicit_rk_butcher<butcher_rk_nssp_21<value_t>>;


template <typename value_t=double>
struct butcher_rk_nssp_32 : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_nssp_32()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/3 , 0 , 0 },{ 0 , 1 , 0 }
    }}, // A
    { 1/2 , 0 , 1/2 }, // b
    { 0 , 1/3 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_nssp_32 = runge_kutta::explicit_rk_butcher<butcher_rk_nssp_32<value_t>>;


template <typename value_t=double>
struct butcher_rk_nssp_33 : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_nssp_33()
  : base_t(
    {{
        { 0 , 0 , 0 },{ -4/9 , 0 , 0 },{ 7/6 , -1/2 , 0 }
    }}, // A
    { 1/4 , 0 , 3/4 }, // b
    { 0 , -4/9 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_nssp_33 = runge_kutta::explicit_rk_butcher<butcher_rk_nssp_33<value_t>>;


template <typename value_t=double>
struct butcher_rk_nssp_53 : public butcher_tableau<5,value_t>
{
  using base_t = butcher_tableau<5,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_nssp_53()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 },{ 1/7 , 0 , 0 , 0 , 0 },{ 0 , 3/16 , 0 , 0 , 0 },{ 0 , 0 , 1/3 , 0 , 0 },{ 0 , 0 , 0 , 2/3 , 0 }
    }}, // A
    { 1/4 , 0 , 0 , 0 , 3/4 }, // b
    { 0 , 1/7 , 3/16 , 1/3 , 2/3 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_nssp_53 = runge_kutta::explicit_rk_butcher<butcher_rk_nssp_53<value_t>>;


template <typename value_t=double>
struct butcher_rk_spp_43 : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_spp_43()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/2 , 0 , 0 , 0 },{ 1/2 , 1/2 , 0 , 0 },{ 1/6 , 1/6 , 1/6 , 0 }
    }}, // A
    { 1/6 , 1/6 , 1/6 , 1/2 }, // b
    { 0 , 1/2 , 1 , 1/2 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_spp_43 = runge_kutta::explicit_rk_butcher<butcher_rk_spp_43<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_22_heun : public butcher_tableau<2,value_t>
{
  using base_t = butcher_tableau<2,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_22_heun()
  : base_t(
    {{
        { 0 , 0 },{ 1 , 0 }
    }}, // A
    { 1/2 , 1/2 }, // b
    { 0 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_22_heun = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_22_heun<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_32 : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_32()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1/2 , 0 , 0 },{ 1/2 , 1/2 , 0 }
    }}, // A
    { 1/3 , 1/3 , 1/3 }, // b
    { 0 , 1/2 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_32 = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_32<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_33 : public butcher_tableau<3,value_t>
{
  using base_t = butcher_tableau<3,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_33()
  : base_t(
    {{
        { 0 , 0 , 0 },{ 1 , 0 , 0 },{ 1/4 , 1/4 , 0 }
    }}, // A
    { 1/6 , 1/6 , 2/3 }, // b
    { 0 , 1 , 1/2 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_33 = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_33<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_42 : public butcher_tableau<4,value_t>
{
  using base_t = butcher_tableau<4,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_42()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 },{ 1/3 , 0 , 0 , 0 },{ 1/3 , 1/3 , 0 , 0 },{ 1/3 , 1/3 , 1/3 , 0 }
    }}, // A
    { 1/4 , 1/4 , 1/4 , 1/4 }, // b
    { 0 , 1/3 , 2/3 , 1 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_42 = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_42<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_53 : public butcher_tableau<5,value_t>
{
  using base_t = butcher_tableau<5,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_53()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 },{ 0.377268915117100 , 0 , 0 , 0 , 0 },{ 0.377268915117100 , 0.377268915117100 , 0 , 0 , 0 },{ 0.163522940897710 , 0.163522940897710 , 0.163522940897710 , 0 , 0 },{ 0.149040593948560 , 0.148312733847240 , 0.148312733847240 , 0.342176968500080 , 0 }
    }}, // A
    { 0.197075963844810 , 0.117803165097650 , 0.117097251937720 , 0.270158749342510 , 0.297864870101040 }, // b
    { 0 , 0.377268915117100 , 0.754537830234190 , 0.490568822693140 , 0.787843030143110 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_53 = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_53<value_t>>;


template <typename value_t=double>
struct butcher_rk_ssp_54 : public butcher_tableau<5,value_t>
{
  using base_t = butcher_tableau<5,value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_rk_ssp_54()
  : base_t(
    {{
        { 0 , 0 , 0 , 0 , 0 },{ 0.391752227003920 , 0 , 0 , 0 , 0 },{ 0.217669096338210 , 0.368410592629590 , 0 , 0 , 0 },{ 0.0826920867095000 , 0.139958502069990 , 0.251891774247380 , 0 , 0 },{ 0.0679662837032000 , 0.115034698444380 , 0.207034898649290 , 0.544974750212370 , 0 }
    }}, // A
    { 0.146811876186610 , 0.248482909245560 , 0.104258830366500 , 0.274438900919600 , 0.226007483193950 }, // b
    { 0 , 0.391752227003920 , 0.586079688967790 , 0.474542363026870 , 0.935010631009240 }  // c
  )
  {}
};
template <typename value_t=double>
using rk_ssp_54 = runge_kutta::explicit_rk_butcher<butcher_rk_ssp_54<value_t>>;



} // namespace ode::butcher