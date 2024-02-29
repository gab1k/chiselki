#ifndef CHEBYSHEV_POLYNOMIALS_HPP_
#define CHEBYSHEV_POLYNOMIALS_HPP_

#include <vector>

namespace adaai {
    // {-1, 0, 3, 0, 0} = 3x^2 - 1
    class ChebyshevPolynomials {
    private:
        std::vector<std::vector<int>> polynoms;

        void count_next() {
            int n = polynoms.size();
            std::vector<int> next(n + 2, 0);
            for (int i = 0; i < n + 1; i++) {
                next[i + 1] = 2 * polynoms[n - 1][i];
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

        explicit ChebyshevPolynomials(int k) {
            polynoms = {{1},
                        {0, 1}};
            for (int i = 0; i < k - 1; i++) {
                count_next();
            }
        }

        std::vector<int> get_polynomial(int k) {
            while (polynoms.size() <= k) {
                count_next();
            }
            return polynoms[k];
        }
    };
}

#endif // CHEBYSHEV_POLYNOMIALS_HPP_