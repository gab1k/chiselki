#ifndef DIFFERENTIATOR_HPP_
#define DIFFERENTIATOR_HPP_

#include <cmath>

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

class AAD22 {
private:
    double m_val; // value of arg
    double m_d1[2]; // df / dx, df/ dy in some point (x, y) in witch m_val computed
    double m_d2[3]; // d^2f/dx^2, d^f/dy^2, d^2f/dxdy
    constexpr AAD22(int i, double v) : m_val(v), m_d2{0, 0, 0},
                                       m_d1{(i == 0 ? 1.0 : 0.0), (i == 0 ? 0.0 : 1.0)} {}; // i == 0 for x, i = 1 for y
public:
    AAD22() = delete;

    constexpr AAD22(double c) : m_val(c), m_d1{0, 0}, m_d2{0, 0, 0} {};

    double get_derivative(WhichD type) {
        if (type == WhichD::x) {
            return m_d1[0];
        } else if (type == WhichD::y) {
            return m_d1[1];
        } else if (type == WhichD::xx) {
            return m_d2[0];
        } else if (type == WhichD::yy) {
            return m_d2[1];
        } else if (type == WhichD::xy) {
            return m_d2[2];
        }
        throw "incorrect type of differentiation";
    }

    constexpr static AAD22 X(double v) {
        return AAD22(0, v);
    }

    constexpr static AAD22 Y(double v) {
        return AAD22(1, v);
    }

    [[nodiscard]] AAD22 my_sin() const {
        AAD22 res = *this;
        res.m_val = sin(this->m_val);
        res.m_d1[0] = this->m_d1[0] * cos(this->m_val);
        res.m_d1[1] = this->m_d1[1] * cos(this->m_val);
        res.m_d2[0] = this->m_d2[0] * cos(this->m_val) - this->m_d1[0] * sin(this->m_val);
        res.m_d2[1] = this->m_d2[1] * cos(this->m_val) - this->m_d1[1] * sin(this->m_val);
        res.m_d2[2] = this->m_d2[2] * cos(this->m_val) - this->m_d1[1] * sin(this->m_val); // CHECK THIS!!!
        return res;
    }

    AAD22 operator+(AAD22 const &r) const {
        AAD22 res = *this;
        res.m_val += r.m_val;
        res.m_d1[0] += r.m_d1[0];
        res.m_d1[1] += r.m_d1[1];
        res.m_d2[0] += r.m_d2[0];
        res.m_d2[1] += r.m_d2[2];
        res.m_d2[2] += r.m_d2[2];
        return res;
    }

    AAD22 operator*(AAD22 const &r) const {
        AAD22 res = *this;
        res.m_val = r.m_val * this->m_val;
        res.m_d1[0] = r.m_d1[0] * this->m_val + r.m_val * this->m_d1[0];
        res.m_d1[1] = r.m_d1[1] * this->m_val + r.m_val * this->m_d1[1];
        res.m_d2[0] = r.m_d2[0] * this->m_val + this->m_d2[0] * r.m_val + 2 * this->m_d1[0] * r.m_d1[0];
        res.m_d2[1] = r.m_d2[1] * this->m_val + this->m_d2[1] * r.m_val + 2 * this->m_d1[1] * r.m_d1[1];;
        res.m_d2[2] = r.m_d2[2] * this->m_val + this->m_d2[2] * r.m_val + r.m_d1[0] * this->m_d1[1] +
                      this->m_d1[0] * r.m_d2[1];
        return res;
    }

    AAD22 operator*(const double &n) const {
        AAD22 res = *this;
        res.m_val *= n;
        res.m_d1[0] *= n;
        res.m_d1[1] *= n;
        res.m_d2[0] *= n;
        res.m_d2[1] *= n;
        res.m_d2[2] *= n;
        return res;
    }


};

AAD22 sin(AAD22 const &val) {
    return val.my_sin();
}

AAD22 operator*(const double &n, const AAD22 &val) {
    return val * n;
}


namespace adaai {
    template<WhichD D, DiffMethod M, typename Callable>
    // Order = 1 or Order = 2;
    double Differentiator(Callable const &F, double x, double y, double h_increase = 1) { // F(x, y) -> a
        if constexpr (M == DiffMethod::FwdAAD) {
            AAD22 X = AAD22::X(x);
            AAD22 Y = AAD22::Y(y);
            AAD22 Res = F(X, Y);
            return Res.get_derivative(D);
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
                return (F(x + h, y + h) - F(x + h, y) - F(x, y + h) + F(x, y)) / (h * h); // CHECK THIS!

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
