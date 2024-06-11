#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "rho.hpp"

class RHS {
public:
    constexpr static int N = 4;

    void operator()(const std::vector<double> &a_y, std::vector<double> &a_rhs, double h = 1e-4) {
        std::vector<double> res = get_diff(a_y);
        for (int i = 0; i < N; i++) {
            a_rhs[i] = res[i];
        }
    }

};

template<typename RHS>
class TimeStepper_RKF45 {
private:
    RHS *m_rhs;
    // нулевой индекс не нужен, нумерация с 1
    const double A[7] = {0, 1.0 / 2, 1.0 / 2, 1, 1, 2.0 / 3, 1.0 / 5};
    const double B[7][7] = {{0, 0,          0,         0,           0,          0,            0},
                            {0, 0,          0,         0,           0,          0,            0},
                            {0, 1.0 / 2,    0,         0,           0,          0,            0},
                            {0, 1.0 / 4,    1.0 / 4,   0,           0,          0,            0},
                            {0, 0,          -1,        2,           0,          0,            0},
                            {0, 7.0 / 27,   10.0 / 27, 0,           1.0 / 27,   0,            0},
                            {0, 28.0 / 625, -1.0 / 5,  546.0 / 625, 54.0 / 625, -378.0 / 625, 0}};
    const double CH[7] = {0, 1.0 / 24, 0, 0, 5.0 / 48, 27.0 / 56, 125.0 / 336};
    const double CT[7] = {0, 1.0 / 8, 0, 2.0 / 3, 1.0 / 16, -27.0 / 56, -125.0 / 336};
public:
    constexpr static int N = RHS::N;

    explicit TimeStepper_RKF45(RHS *a_rhs) : m_rhs(a_rhs) {}

    std::pair<double, double>
    operator()(double a_t, const std::vector<double> &a_y, double h,
               std::vector<double> &a_y_next, double eps = 1e-6) const { // return next_t, next_h
        // RKF45 for 1 step
        std::vector<double> next_diff(N);

        // fill k1
        std::vector<double> k1(N);
        (*m_rhs)(a_y, k1, a_t + A[1] * h);
        for (int i = 0; i < N; i++) {
            k1[i] *= h;
            next_diff[i] = a_y[i] + B[2][1] * k1[i];  // point to find k2
        }

        // fill k2
        std::vector<double> k2(N);
        (*m_rhs)(next_diff, k2, a_t + A[2] * h);
        for (int i = 0; i < N; i++) {
            k2[i] *= h;
            next_diff[i] = a_y[i] + B[3][1] * k1[i] + B[3][2] * k2[i]; // point to find k3
        }

        // fill k3
        std::vector<double> k3(N);
        (*m_rhs)(next_diff, k3, a_t + A[3] * h);
        for (int i = 0; i < N; i++) {
            k3[i] *= h;
            next_diff[i] = a_y[i] + B[4][1] * k1[i] + B[4][2] * k2[i] + B[4][3] * k3[i]; // point to find k4
        }

        // fill k4
        std::vector<double> k4(N);
        (*m_rhs)(next_diff, k4, a_t + A[4] * h);
        for (int i = 0; i < N; i++) {
            k4[i] *= h;
            next_diff[i] =
                    a_y[i] + B[5][1] * k1[i] + B[5][2] * k2[i] + B[5][3] * k3[i] + B[5][4] * k4[i]; // point to find k5
        }

        // fill k5
        std::vector<double> k5(N);
        (*m_rhs)(next_diff, k5, a_t + A[5] * h);
        for (int i = 0; i < N; i++) {
            k5[i] *= h;
            next_diff[i] =
                    a_y[i] + B[6][1] * k1[i] + B[6][2] * k2[i] + B[6][3] * k3[i] + B[6][4] * k4[i] + B[6][5] * k5[i];
        }

        // fill k6
        std::vector<double> k6(N);
        (*m_rhs)(next_diff, k6, a_t * A[6] * h);
        for (double &elem: k6) {
            elem *= h;
        }

        // find step value
        for (int i = 0; i < N; i++) {
            a_y_next[i] =

                    a_y[i] + CH[1] * k1[i] + CH[2] * k2[i] + CH[3] * k3[i] + CH[4] * k4[i] + CH[5] * k5[i] +
                    CH[6] * k6[i];
        }

        // find next_h
        double max_err = 1e-10;
        for (int i = 0; i < N; i++) {
            max_err = fmax(max_err, fabs(CT[1] * k1[i] + CT[2] * k2[i] + CT[3] * k3[i] + CT[4] * k4[i] + CT[5] * k5[i] +
                                         CT[6] * k6[i]));
        }
        double h_new = 0.9 * h * pow(eps / max_err, 1.0 / 5);
        return {a_t + h_new, h_new};
    }
};


template<typename Stepper, typename Observer>
class ODE_Integrator {
private:
    Stepper const *const m_stepper;
    Observer *const m_observer;
public:
    void operator()(double a_t0, const std::vector<double> &a_y0, double a_tEnd, std::vector<double> &a_yEnd) {
//        static_assert(a_y0.size() == Stepper::N);
//        static_assert(a_yEnd.size() == Stepper::N);
        double curr_t = a_t0;
        double curr_h = 1e-5;
        std::vector<double> curr_y = a_y0;
        for (int i = 1; curr_t < a_tEnd; i++) {
            std::pair<double, double> res = (*m_stepper)(curr_t, curr_y, curr_h, a_yEnd);
            curr_t = res.first;
            curr_h = res.second;
            for (int j = 0; j < a_yEnd.size(); j++) {
                curr_y[j] = a_yEnd[j];
            }
            if (!(*m_observer)(curr_t, i, curr_y)) {
                break;
            }
        }
    }

    ODE_Integrator(Stepper *a_stepper, Observer *a_observer) : m_stepper(a_stepper), m_observer(a_observer) {}
};


class ObserverMaxRange {
private:
    double max_x = 0;
public:
    bool operator()(double a_curr_t, int a_n, std::vector<double> &a_curr_y) {
//        std::cout << "Step number = " << a_n << ". Current time = " << a_curr_t << ". Current x = " << a_curr_y[0]
//                  << ". Current y = " << a_curr_y[2]
//                  << "\n";
        if (a_curr_y[2] <= 0 && a_n > 1)return false;
        max_x = fmax(max_x, a_curr_y[0]);
        return true;
    }

    [[nodiscard]] double get_max_distance() const {
        return max_x;
    }

    void set_zero_distance() {
        max_x = 0;
    }
};

class ObserverBestCorner {
private:
    std::vector<std::vector<double>> points; // t, x, y
public:
    ObserverBestCorner() : points(3, std::vector<double>()) {}

    bool operator()(double a_curr_t, int a_n, std::vector<double> &a_curr_y) {
//        std::cout << "Step number = " << a_n << ". Current time = " << a_curr_t << ". Current x = " << a_curr_y[0]
//                  << ". Current y = " << a_curr_y[2]
//                  << "\n";
        if (a_curr_y[2] <= 0 && a_n > 1)return false;
        points[0].push_back(a_curr_t); // t
        points[1].push_back(a_curr_y[0]); // x
        points[2].push_back(a_curr_y[2]); // y
        return true;
    }

    void fill_point(std::vector<double> &t, std::vector<double> &x, std::vector<double> &y) {
        t.clear(), x.clear(), y.clear();
        for (int i = 0; i < points[0].size(); i++) {
            t.push_back(points[0][i]);
            x.push_back(points[1][i]);
            y.push_back(points[2][i]);
        }
    }
};

void print_points_to_plot(std::vector<double> &x, std::vector<double> &y) {
    std::string filePath = "../atmosphere_model/logs/points_logs.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out);

    ofs << "x = [\n";
    for (double i: x) {
        ofs << i << ", ";
    }
    ofs << "\n]\ny = [\n";
    for (double i: y) {
        ofs << i << ", ";
    }
    ofs << "\n]\n";
    ofs.close();
}

int main() {
    RHS r;
    TimeStepper_RKF45<RHS> stepper(&r);
    ObserverMaxRange maxRange;
    ODE_Integrator integrator(&stepper, &maxRange);
    double v0 = 1640;
    std::pair<double, double> res = {0, 0}; // alpha, max_height
    for (double alpha = 40.0; alpha <= 60; alpha += 1) {
        double radian = alpha * M_PI / 180.0;
        std::vector<double> end(RHS::N), start = {0, v0 * cos(radian), 0, v0 * sin(radian)}; // x, v_x, y, v_y

        maxRange.set_zero_distance();
        integrator(0, start, 1e9, end);
        if (res.second < maxRange.get_max_distance()) {
            res = {alpha, maxRange.get_max_distance()};
        }
    }
    std::cout << "max distance = " << res.second << ". With corner = " << res.first << "\n";
    // best corner = 53

    ObserverBestCorner get_plot;
    ODE_Integrator best_corner_integrator(&stepper, &get_plot);
    double alpha = 53 * M_PI / 180;
    std::vector<double> t, x, y, end(RHS::N), start = {0, v0 * cos(alpha), 0, v0 * sin(alpha)};
    best_corner_integrator(0, start, 1e9, end);
    get_plot.fill_point(t, x, y);
    print_points_to_plot(x, y);
    return 0;
}
