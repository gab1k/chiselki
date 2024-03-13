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

#endif // ENUMS_HPP_