#ifndef ENUMS_HPP_
#define ENUMS_HPP_

#include <unordered_map>
#include <string>

enum class DiffMethod : int {
    stencil3,
    stencil3Extra,
    stencil5,
    stencil5Extra,
    FwdAAD
};

enum class WhichD : int {
    x, // df / dx
    y, // df / dy
    xx, // df / dxdx
    yy, // df / dydy
    xy // df / dxdy
};

std::unordered_map<WhichD, std::string> namesT{
        {WhichD::x,  "X"},
        {WhichD::y,  "Y"},
        {WhichD::xx, "XX"},
        {WhichD::yy, "YY"},
        {WhichD::xy, "XY"}
};
std::unordered_map<DiffMethod, std::string> namesM{
        {DiffMethod::stencil3,      "stencil3"},
        {DiffMethod::stencil3Extra, "stencil3Extra"},
        {DiffMethod::stencil5,      "stencil5"},
        {DiffMethod::stencil5Extra, "stencil5Extra"},
        {DiffMethod::FwdAAD,        "FwdAAD"}
};

#endif // ENUMS_HPP_