#ifndef DIFFERENTIATOR_HPP_
#define DIFFERENTIATOR_HPP_

#include <cmath>
#include "aad22.hpp"
#include "enums.hpp"


namespace adaai {
    template<WhichD D, DiffMethod M, typename Callable>
    double Differentiator(Callable const &F, double x, double y, double h_increase = 1) { // F(x, y) -> a
        if constexpr (M == DiffMethod::FwdAAD) {
            try {
                AAD22 X = AAD22::X(x);
                AAD22 Y = AAD22::Y(y);
                AAD22 Res = F(X, Y);
                return Res.get_derivative(D);
            } catch (...) {
                throw "Callable function is not AAD22 class";
            }
        } else {
            double h = 1e-4 * h_increase;
            if (std::abs(x) >= 1) {
                h *= std::abs(x);
            }
            if constexpr (M == DiffMethod::stencil3) {
                if constexpr (D == WhichD::x) {
                    return (F(x + h, y) - F(x - h, y)) / (2 * h);
                } else if constexpr (D == WhichD::y) {
                    return (F(x, y + h) - F(x, y - h)) / (2 * h);
                } else if constexpr (D == WhichD::xx) {
                    return (F(x + h, y) - 2 * F(x, y) + F(x - h, y)) / (h * h);
                } else if constexpr (D == WhichD::yy) {
                    return (F(x, y + h) - 2 * F(x, y) + F(x, y - h)) / (h * h);
                }
                double d_p = (F(x + h, y + h) - F(x - h, y + h)) / (2 * h);
                double d_m = (F(x + h, y - h) - F(x - h, y - h)) / (2 * h);
                return (d_p - d_m) / (2 * h);

            } else if constexpr (M == DiffMethod::stencil5) {
                if constexpr (D == WhichD::x) {
                    return (-F(x + 2 * h, y) + 8 * F(x + h, y) - 8 * F(x - h, y) + F(x - 2 * h, y)) / (12 * h);
                } else if constexpr (D == WhichD::y) {
                    return (-F(x, y + 2 * h) + 8 * F(x, y + h) - 8 * F(x, y - h) + F(x, y - 2 * h)) / (12 * h);
                } else if constexpr (D == WhichD::xx) {
                    return (-F(x - 2 * h, y) + 16 * F(x - h, y) - 30 * F(x, y) + 16 * F(x + h, y) - F(x + 2 * h, y)) /
                           (12 * h * h);
                } else if constexpr (D == WhichD::yy) {
                    return (-F(x, y - 2 * h) + 16 * F(x, y - h) - 30 * F(x, y) + 16 * F(x, y + h) - F(x, y + 2 * h)) /
                           (12 * h * h);
                }
                double d_p = (F(x + h, y + h) - F(x - h, y + h)) / (2 * h);
                double d_m = (F(x + h, y - h) - F(x - h, y - h)) / (2 * h);
                return (d_p - d_m) / (2 * h);

            } else if constexpr (M == DiffMethod::stencil3Extra) {
                double d3_h = Differentiator<D, DiffMethod::stencil3>(F, x, y);
                double d3_h2 = Differentiator<D, DiffMethod::stencil3>(F, x, y, 0.5); // n = 2;
                return (4 * d3_h2 - d3_h) / 3;
            } else if constexpr (M == DiffMethod::stencil5Extra) {
                double d3_h = Differentiator<D, DiffMethod::stencil5>(F, x, y);
                double d3_h2 = Differentiator<D, DiffMethod::stencil5>(F, x, y, 0.5); // n = 2;
                return (4 * d3_h2 - d3_h) / 3;
            }
        }

    }

} // namespace adaai

#endif // DIFFERENTIATOR_HPP_
