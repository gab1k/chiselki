#include <vector>
#include <iostream>
#include <math.h>

class RHS {
public:
    void operator()(const std::vector<double> &a_y, std::vector<double> &a_rhs, double _t0, double _t) {
        // and solve from a_y, to a_rhs. Change values of a_rhs
        // a_y -> R ^ {2N + 1}
        // a_rhs -> R ^ {N}
    }
};


template<typename RHS>
class TimeStepper_everhart {
public:
    constexpr static int N = RHS::N;
    constexpr static int K = 5;

    TimeStepper_everhart(RHS *a_rhs) : rhs(a_rhs) {
        t.resize(K + 1);
        F.resize(K + 1, std::vector<std::vector<double>>(K + 1, std::vector<double>(N)));
    }

    std::pair<double, double> operator()(double a_t, const std::vector<double> &a_y, double h,
                                         std::vector<double> &a_y_next, double eps = 1e-6) {
        rhs(a_y, F0, a_t, a_t); // find f(t_0, y(t_0), dy/dt(t_0);
        find_divided_difference(a_y, a_t);
        // a_y.size() = N * 2 + 1
        // a_y[0] = t
        // a_y[1, ..., N] = y
        // a_y[N+1, ..., 2N] = dy/dt
        for (int i = 0; i <= K; i++) {
            t[i] = a_t + h * i / K;
        }
    }

private:
    RHS *rhs;
    std::vector<double> t;
    std::vector<std::vector<std::vector<double>>> F; // index start, size, vector F
    // if F[2][3] its means F[t2, t3, t4]. F[1][1] its means F[1]
    std::vector<double> F0;
//    std::vector<double> start, res;

    void find_divided_difference(const std::vector<double> &a_y, double _t0) { //
        std::vector<double> F_i(N);
        std::vector<double> curr_a_y(2 * N + 1);
        for (int i = 0; i <= K; i++) {
            curr_a_y[0] = t[i];
            fill_dy_dt_for_divided_diff(a_y, t[i], _t0, curr_a_y);
            fill_y_for_divided_diff(a_y, t[i], _t0, curr_a_y);

            rhs(curr_a_y, F_i, _t0, t[i]);
            F[i][1] = F_i;
        }
        for (int i = 0; i < K; i++) {
            for (int j = 2; i + j - 1 <= K; j++) {
                // F[i][j] = (F[i+1][j-1] - F[i][j-1]) / (t[i+j-1] - t[i])
                for (int id; id < N; id++) {
                    F[i][j][id] = (F[i + 1][j - 1][id] - F[i][j - 1][id]) / (t[i + j - 1] - t[i]);
                }
            }
        }
    }

    void fill_dy_dt_for_divided_diff(const std::vector<double> &a_y, double _t, double _t0,
                                     std::vector<double> &where_write) {
        std::vector<double> dy_dt(N);
        for (int i = 0; i < N; i++) {
            where_write[N + 1 + i] = a_y[N + i + 1] + F0[i] * (_t - _t0); // dy/dt = dy/dt + F0 * (t-t0)
        }
    }

    std::vector<double>
    fill_y_for_divided_diff(const std::vector<double> &a_y, double _t, double _t0, std::vector<double> &where_write) {
        std::vector<double> y(N);
        for (int i = 0; i < N; i++) {
            where_write[i + 1] = a_y[1 + i] + (_t - _t0) * a_y[N + i + 1] + F0[i] * pow((_t - _t0), 2) / 2;
        }
    }
};
