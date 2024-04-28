#include <vector>
#include <iostream>
#include <cmath>

class RHS {
public:
    constexpr static int N = 3;

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
    constexpr static int COUNT_APPROACH = 10;

    TimeStepper_everhart(RHS *a_rhs) : rhs(a_rhs) {
        t.resize(K + 1);
        B.resize(K + 1);
        F.resize(K + 1, std::vector<std::vector<double>>(K + 1, std::vector<double>(N)));
    }

    std::pair<double, double> operator()(double a_t, const std::vector<double> &a_y, double h,
                                         std::vector<double> &a_y_next, double eps = 1e-6) { // return next_t, next_h
        // a_y.size() = N * 2 + 1
        // a_y[0] = t
        // a_y[1, ..., N] = y
        // a_y[N+1, ..., 2N] = dy/dt

        // a_y_next.size() = N = (dy)^2 / (dt)^2

        for (int i = 0; i <= K; i++) {
            t[i] = a_t + h * i / K;
        }
        rhs(a_y, F0, a_t, a_t); // find f(t_0, y(t_0), dy/dt(t_0);
        for (int step = 0; step < COUNT_APPROACH; step++) {
            find_divided_difference(a_y, a_t, step == 0);
            find_B();
        }
        std::vector<double> rhs_val(2 * N + 1); // t, y, dy/dt
        rhs_val[0] = a_t + h;
        fill_y(a_y, a_t + h, a_t, rhs_val);
        fill_dy_dt(a_y, a_t + h, a_t, rhs_val);
        rhs(rhs_val, a_y_next, a_t, a_t + h);
        return {a_t + h, h};
    }

private:
    RHS *rhs;
    std::vector<double> t;
    std::vector<std::vector<std::vector<double>>> F; // index start, size, vector F
    // if F[2][3] its means F[t2, t3, t4]. F[1][1] its means F[1]
    std::vector<double> F0, B;
//    std::vector<double> start, res;

    void find_B() {
/*
 L(t) = F(t0) + \sum\limits_{i=1}^k * (t - t0) * ... * (t - t_{i-1})
 и мы знаем, что k = 5, t_i = t0 + h/k * i = t0 + ih/5
 То, чего не было на практике, но чем мы пользуемся:
 L(t) =  F[t0] +
         F[t0, t1] * (t - t0) +
         F[t0, t1, t2] * (t - t0) * (t - t1) +
         F[t0, t1, t2, t3] * (t - t0) * (t - t1) * (t - t2) +
         F[t0, t1, t2, t3, t4) *  (t - t0) * (t - t1) * (t - t2) * (t - t3) +
         F[t0, t1, t2, t3, t4, t5) *  (t - t0) * (t - t1) * (t - t2) * (t - t3) * (t - t4)
     = (подставим t_i = t0 + ih/5) =
         F[t0) +
         F[t0, t1] * (t - t0) +
         F[t0, t1, t2] * (t^2 - 2*t0*t + t0^2 + ht/5 - h*t0/5) +
         и так далее, перписывать не стану, просто соберу по стеменям t:
    = Пусть a = F[t0], b = F[t0, t1], ..., f = F[t0, ..., t5] =
        Тогда L(t) = \sum\limits_{j=0}^k (B[j] * (t - t0)^j)
        B[j] - выражается через a,b,c,d,e,f,h,t0
 */

    }

    void find_divided_difference(const std::vector<double> &a_y, double _t0, bool initial = false) { //
        std::vector<double> F_i(N);
        std::vector<double> curr_a_y(2 * N + 1);
        for (int i = 0; i <= K; i++) {
            curr_a_y[0] = t[i];
            fill_dy_dt(a_y, t[i], _t0, curr_a_y, initial);
            fill_y(a_y, t[i], _t0, curr_a_y, initial);

            rhs(curr_a_y, F_i, _t0, t[i]);
            F[i][1] = F_i;
        }
        for (int i = 0; i < K; i++) {
            for (int j = 2; i + j - 1 <= K; j++) {
                // F[i][j] = (F[i+1][j-1] - F[i][j-1]) / (t[i+j-1] - t[i]) (вообще можно сказать что (t[i+j-1] - t[i]) = k * (j - 1) / 5)
                for (int id; id < N; id++) {
                    F[i][j][id] = (F[i + 1][j - 1][id] - F[i][j - 1][id]) / (t[i + j - 1] - t[i]);
                }
            }
        }
    }

    void fill_dy_dt(const std::vector<double> &a_y, double _t, double _t0,
                    std::vector<double> &where_write, bool initial = false) {
        std::vector<double> dy_dt(N);
        if (initial) {
            for (int i = 0; i < N; i++) {
                where_write[N + 1 + i] = a_y[N + i + 1] + F0[i] * (_t - _t0); // dy/dt = dy/dt + F0 * (t-t0)
            }
        } else {
            for (int i = 0; i < N; i++) {
                where_write[N + 1 + i] = a_y[N + i + 1];
                for (int j = 0; j <= K; j++) {
                    where_write[N + 1 + i] += (B[j] * pow(_t - _t0, j + 1)) / (j + 1);
                }
            }
        }
    }

    std::vector<double>
    fill_y(const std::vector<double> &a_y, double _t, double _t0, std::vector<double> &where_write,
           bool initial = false) {
        std::vector<double> y(N);
        if (initial) {
            for (int i = 0; i < N; i++) {
                where_write[i + 1] = a_y[1 + i] + (_t - _t0) * a_y[N + i + 1] + F0[i] * pow((_t - _t0), 2) / 2;
            }
        } else {
            for (int i = 0; i < N; i++) {
                where_write[i + 1] = a_y[i + 1] + a_y[N + i + 1] * (_t - _t0);
                for (int j = 0; j <= K; j++) {
                    where_write[i + 1] += (B[j] * pow(_t - _t0, j + 2)) / ((j + 1) * (j + 2));
                }
            }
        }
    }
};

int main() {
    RHS rhs;
    TimeStepper_everhart stepper(&rhs);
}