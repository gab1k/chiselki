#ifndef AAD22_HPP_
#define AAD22_HPP_

#include "enums.hpp"

class AAD22 {
private:
    double m_val;
    double m_d1[2];
    double m_d2[3];
    constexpr AAD22(int i, double v) : m_val(v), m_d2{0, 0, 0},
                                       m_d1{(i == 0 ? 1.0 : 0.0), (i == 0 ? 0.0 : 1.0)} {}; // i == 0 for x, i = 1 for y
public:
    AAD22() = delete;

    constexpr explicit AAD22(double c) : m_val(c), m_d1{0, 0}, m_d2{0, 0, 0} {};

    constexpr static AAD22 X(double v) {
        return {0, v};
    }

    constexpr static AAD22 Y(double v) {
        return {1, v};
    }

    double get_derivative(WhichD type);

    [[nodiscard]] AAD22 my_sin() const;

    AAD22 operator+(AAD22 const &r) const;

    AAD22 operator+(double const &c) const;


    AAD22 operator*(AAD22 const &r) const;

    AAD22 operator*(double const&n) const;


};

AAD22 operator*(double const&n, AAD22 const&val);

AAD22 operator+(double const&n, AAD22 const&val);


AAD22 sin(AAD22 const &val);

#endif // AAD22_HPP_

