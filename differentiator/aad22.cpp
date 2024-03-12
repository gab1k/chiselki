#include "aad22.hpp"
#include <cmath>

AAD22 AAD22::my_sin() const {
    AAD22 res = *this;
    res.m_val = sin(this->m_val);
    res.m_d1[0] = this->m_d1[0] * cos(this->m_val);
    res.m_d1[1] = this->m_d1[1] * cos(this->m_val);
    res.m_d2[0] = this->m_d2[0] * cos(this->m_val) - this->m_d1[0] * this->m_d1[0] * sin(this->m_val);
    res.m_d2[1] = this->m_d2[1] * cos(this->m_val) - this->m_d1[1] * this->m_d1[1] * sin(this->m_val);
    res.m_d2[2] = this->m_d2[2] * cos(this->m_val) - this->m_d1[1] * this->m_d1[0] * sin(this->m_val);
    return res;
}

AAD22 sin(AAD22 const &val) {
    return val.my_sin();
}


AAD22 AAD22::my_cos() const {
    AAD22 res = *this;
    res.m_val = cos(this->m_val);
    res.m_d1[0] = -this->m_d1[0] * sin(this->m_val);
    res.m_d1[1] = -this->m_d1[1] * sin(this->m_val);
    res.m_d2[0] = -this->m_d1[0] * this->m_d1[0] * cos(this->m_val) - this->m_d2[0] * sin(this->m_val);
    res.m_d2[1] = -this->m_d1[1] * this->m_d1[1] * cos(this->m_val) - this->m_d2[1] * sin(this->m_val);
    res.m_d2[2] = -this->m_d1[0] * this->m_d1[1] * cos(this->m_val) - this->m_d2[2] * sin(this->m_val);
    return res;
}

AAD22 cos(AAD22 const &val) {
    return val.my_cos();
}

AAD22 AAD22::my_exp() const {
    AAD22 res = *this;
    res.m_val = exp(this->m_val);
    res.m_d1[0] = exp(this->m_val) * this->m_d1[0];
    res.m_d1[1] = exp(this->m_val) * this->m_d1[1];
    res.m_d2[0] = exp(this->m_val) * (this->m_d1[0] * this->m_d1[0] + this->m_d2[0]);
    res.m_d2[1] = exp(this->m_val) * (this->m_d1[1] * this->m_d1[1] + this->m_d2[1]);
    res.m_d2[2] = exp(this->m_val) * (this->m_d1[0] * this->m_d1[1] + this->m_d2[2]);
    return res;
}

AAD22 exp(AAD22 const &val) {
    return val.my_exp();
}


AAD22 AAD22::operator+(AAD22 const &r) const {
    AAD22 res = *this;
    res.m_val += r.m_val;
    res.m_d1[0] += r.m_d1[0];
    res.m_d1[1] += r.m_d1[1];
    res.m_d2[0] += r.m_d2[0];
    res.m_d2[1] += r.m_d2[1];
    res.m_d2[2] += r.m_d2[2];
    return res;
}

AAD22 AAD22::operator+(const double &c) const {
    return *this + AAD22(c);
}

AAD22 operator+(double const &n, AAD22 const &val) {
    return val + n;
}

AAD22 AAD22::operator+=(const AAD22 &r) {
    *this = *this + r;
    return *this;
}

AAD22 AAD22::operator+=(const double &c) {
    return *this += AAD22(c);
}


AAD22 AAD22::operator*(AAD22 const &r) const {
    AAD22 res = *this;
    res.m_val = r.m_val * this->m_val;
    res.m_d1[0] = r.m_d1[0] * this->m_val + r.m_val * this->m_d1[0];
    res.m_d1[1] = r.m_d1[1] * this->m_val + r.m_val * this->m_d1[1];
    res.m_d2[0] = r.m_d2[0] * this->m_val + this->m_d2[0] * r.m_val + 2 * this->m_d1[0] * r.m_d1[0];
    res.m_d2[1] = r.m_d2[1] * this->m_val + this->m_d2[1] * r.m_val + 2 * this->m_d1[1] * r.m_d1[1];
    res.m_d2[2] = r.m_d2[2] * this->m_val + this->m_d2[2] * r.m_val + r.m_d1[0] * this->m_d1[1] +
                  this->m_d1[0] * r.m_d1[1];
    return res;
}

AAD22 AAD22::operator*(double const &n) const {
    return *this * AAD22(n);
}

AAD22 operator*(double const &n, AAD22 const &val) {
    return val * n;
}


AAD22 AAD22::operator*=(const AAD22 &r) {
    *this = *this * r;
    return *this;
}

AAD22 AAD22::operator*=(const double &c) {
    return *this *= AAD22(c);
}

AAD22 AAD22::operator-(const AAD22 &r) const {
    return *this + (-1 * AAD22(r));
}

AAD22 AAD22::operator-(const double &c) const {
    return *this - AAD22(c);
}

AAD22 operator-(double const &n, AAD22 const &val) {
    return AAD22(n) - val;
}

AAD22 AAD22::operator-=(const AAD22 &r) {
    *this = *this - r;
    return *this;
}

AAD22 AAD22::operator-=(const double &c) {
    return *this -= AAD22(c);
}

AAD22 AAD22::operator/(const AAD22 &r) const {
    AAD22 res = *this;
    res.m_val = this->m_val / r.m_val;
    res.m_d1[0] = (this->m_d1[0] * r.m_val - this->m_val * r.m_d1[0]) / (r.m_val * r.m_val);
    res.m_d1[1] = (this->m_d1[1] * r.m_val - this->m_val * r.m_d1[1]) / (r.m_val * r.m_val);
    res.m_d2[0] =
            (-r.m_val * (2 * this->m_d1[0] * r.m_d1[0] + this->m_val * r.m_d2[0]) + this->m_d2[0] * r.m_val * r.m_val +
             2 * this->m_val * r.m_d1[0] * r.m_d1[0]) / (r.m_val * r.m_val * r.m_val);
    res.m_d2[1] =
            (-r.m_val * (2 * this->m_d1[1] * r.m_d1[1] + this->m_val * r.m_d2[1]) + this->m_d2[1] * r.m_val * r.m_val +
             2 * this->m_val * r.m_d1[1] * r.m_d1[1]) / (r.m_val * r.m_val * r.m_val);
    res.m_d2[2] = (-r.m_val * (this->m_d1[0] * r.m_d1[1] + this->m_d1[1] * r.m_d1[0] + this->m_val * r.m_d2[2]) +
                   this->m_d2[2] * r.m_val * r.m_val + 2 * this->m_val * r.m_d1[0] * r.m_d1[1]) /
                  (r.m_val * r.m_val * r.m_val);
    return res;
}

AAD22 AAD22::operator/(const double &n) const {
    return *this / AAD22(n);
}

AAD22 operator/(double const &n, AAD22 const &val) {
    return AAD22(n) / val;
}

AAD22 AAD22::operator/=(const AAD22 &r) {
    *this = *this / r;
    return *this;
}

AAD22 AAD22::operator/=(const double &c) {
    return *this /= AAD22(c);
}

double AAD22::get_derivative(WhichD type) {
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
