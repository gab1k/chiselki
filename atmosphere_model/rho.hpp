#ifndef GET_RHO
#define GET_RHO

#include <cmath>
#include <vector>

std::vector<double> get_diff(const std::vector<double> &u);

class Rho {
private:
    const double h_l[5], r_l[4], T_l[4], p_l[4];
    const double g = 9.80655;
    const double R_air = 287.0528;

    int get_l(double h) {
        for (int i = 1; i < 5; i++) {
            if (h_l[i] > h) {
                return i - 1;
            }
        }
        // Going beyond the definition
        return 4;
    }

    double get_T_new_layer(int l) {
        if (l == 0)return T_l[0];
        double p_h = get_p_new_layer(l);
        double rho_h = p_h / (R_air * (T_l[l - 1] - r_l[l - 1] * (h_l[l] - h_l[l - 1])));
        return p_h / (rho_h * R_air);
    }

    double get_p_new_layer(int l) {
        if (l == 0) return p_l[0];
        return get_p(h_l[l], -1);
    }

public:
    Rho() : h_l{0, 11000, 20000, 32000, 47000},
            r_l{6.5 * 1e-3, 0, -1e-3, -2.8 * 1e-3},
            T_l{288.15, 216.65, 216.65, 228.65},
            p_l{101325, 100343, 24274.1, 24318.8} {}

    void print_all_coeff() {
        for (int l = 0; l < 4; l++) {
            std::cout << "Layer: " << l << "\n";
            std::cout << "h = " << h_l[l] << ", r = " << r_l[l] << ", " << "T = " << get_T_new_layer(l) << ", "
                      << "p = " << get_p_new_layer(l) << "\n";
        }
    }

    double get_p(double h, int layer_delta = 0) {
        int l = get_l(h) + layer_delta;
        if (r_l[l] == 0) {
            return p_l[l] * std::exp(-(g * (h - h_l[l])) / (R_air * T_l[l]));
        }
        return p_l[l] * std::exp(g / R_air * log(1 - (r_l[l] * (h - h_l[l])) / T_l[l]));
    }

    double operator()(double h) {
        int l = get_l(h);
        return get_p(h) / (R_air * (T_l[l] - r_l[l] * (h - h_l[l])));
    }
};

class CD {
private:
    double Mp = 0.5 / 118.0; // 0.5 M = 118 pixels
    double CDp = 0.1 / 80.0; // 0.1 CD = 80 pixels
    double a = 0.46431097913627956, b = 0.67929407312546719;
    std::vector<double> coef; // coef[i] = (point[i + 1].s - point[i].s) / (point[i + 1].f - point[i].f)
    std::vector<std::pair<double, double>> points{{0.25,           0.1 + CDp * 11},     // M, CD(M)
                                                  {0.5 + Mp * 23,  0.1 + CDp * 14},
                                                  {0.5 + Mp * 70,  0.1 + CDp * 21},
                                                  {0.5 + Mp * 83,  0.1 + CDp * 29},
                                                  {0.5 + Mp * 95,  0.1 + CDp * 55},
                                                  {0.5 + Mp * 106, 0.2 + CDp * 39},
                                                  {1.0,            0.4 + CDp * 32},
                                                  {1.0 + Mp * 12,  0.4 + CDp * 24},
                                                  {1.0 + Mp * 25,  0.4 + CDp * 8},
                                                  {1.0 + Mp * 37,  0.4 + CDp * 4},
                                                  {1.0 + Mp * 49,  0.3 + CDp * 79},
                                                  {1.0 + Mp * 95,  0.3 + CDp * 54},
                                                  {1.5 + Mp * 25,  0.3 + CDp * 29},
                                                  {1.5 + Mp * 70,  0.3 + CDp * 13},
                                                  {2 + Mp * 25,    0.2 + CDp * 61},
                                                  {2 + Mp * 47,    0.2 + CDp * 52}
    };

    int get_i(double M) {
        for (int i = 0; i < points.size(); i++) {
            if (M <= points[i].first) {
                return i - 1;
            }
        }
        return (int) points.size();
    }

public:

    CD() {
        coef.resize(points.size() - 1);
        for (int i = 1; i < points.size(); i++) {
            coef[i - 1] = (points[i].second - points[i - 1].second) / (points[i].first - points[i - 1].first);
        }
    }

    double operator()(double M) {
        int i = get_i(M);
        if (i < points.size() - 1) { // if between 2 points
            return points[i].second + coef[i] * (M - points[i].first);
        }
        return a * pow(M, -b);
    }
};

void find_a_b();

#endif // GET_RHO
