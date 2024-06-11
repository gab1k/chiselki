#ifndef CHEBYSHEV_POLYNOMIALS_HPP_
#define CHEBYSHEV_POLYNOMIALS_HPP_

#include <vector>

namespace adaai {

    // {-1, 0, 3, 0, 0} = 3x^2 - 1
    class ChebyshevPolynomials {
    private:
        std::vector<std::vector<long long>> polynoms;

        void count_next() {
            unsigned long long n = polynoms.size();
            std::vector<long long> next(n + 2, 0);
            for (int i = 0; i < n + 1; i++) {
                next[i + 1] = (polynoms[n - 1][i] << 1);
            }
            for (int i = 0; i < n; i++) {
                next[i] -= polynoms[n - 2][i];
            }
            polynoms.push_back(next);
        }

    public:
        ChebyshevPolynomials() {
            polynoms = {{1},
                        {0, 1}};
        }

        explicit ChebyshevPolynomials(unsigned long long k) {
            polynoms = {{1},
                        {0, 1}};
            for (int i = 0; i < k - 1; i++) {
                count_next();
            }
        }

        std::vector<long long> get_polynomial(unsigned long long k) {
            while (polynoms.size() <= k) {
                count_next();
            }
            return polynoms[k];
        }
    };

} // namespace adaai

#endif // CHEBYSHEV_POLYNOMIALS_HPP_
