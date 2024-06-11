#ifndef EXP_PADE_HPP_
#define EXP_PADE_HPP_

#include <vector>

namespace adaai {

    template<typename F>
    constexpr F exp_pade(F x) {
//        std::vector<F> members = {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0}; // for x^0, x^1, ..., x^5.
//        std::vector<F> members = {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0}; // to x^7
        std::vector<F> members = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0,
                                  110880.0, 3960.0, 90.0, 1.0}; // to x^9

        F x_pow = 1;
        F divisible = 0;
        F divisor = 0;
        for (unsigned pow = 0; pow < members.size(); ++pow) {
            F number = x_pow * members[pow];
            divisible += number;
            if (pow & 1) {
                divisor -= number;
            } else {
                divisor += number;
            }
            x_pow *= x;
        }

        return divisible / divisor;
    }

} // namespace adaai

#endif // EXP_PADE_HPP_
