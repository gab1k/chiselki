#include "stepper.hpp"
#include <string>
#include <fstream>

const double MU = 398600.4; // km^3/sec^2
const double Re = 6378.137; // km радиус Земли
const double J2 = 1.0827 * 1e-3;

class RHS_U2 {
public:
    // r = (x, y, z); |r| = sqrt(x^2+y^2+z^2)

    // sin (phi) = z / |r|
    // U_2(r) = mu/|r| * (1 - J_2 * (R_e/|r|)^2 * (3 * (z/|r|)^2 - 1) / 2 = mu/|r| - (J_2 * mu * R_e^2 * 3 * z^2) / (2 * |r|^5) - (mu * R_e^2) / (2|r|^3)

    // d/dx mu/|r| = (-mu * x) / pow(x^2+y^2+z^2, 3/2)
    // d/dx (J_2 * mu * R_e^2 * 3 * z^2) / (2 * |r|^5) = J_2 * mu * R_e^2 * 3 * z^2 / 2 * (- (5x) / pow(x^2+y^2+z^2, 7/2) )

    // d/dz (J_2 * mu * R_e^2 * 3 * z^2) / (2 * |r|^5) = J_2 * mu * R_e^2 * 3 / 2 * ( 2x^2 * z + 2y^2 * z - 3z^3 / pow(x^2+y^2+z^2, 7/2)

    // d/dx = (mu * R_e^2) / (2|r|^3) = mu * R_e^2 / 2 * (- (3x) /  pow(x^2+y^2+z^2, 5/2))

    constexpr static int N = 3;

    void operator()(const std::vector<double> &a_y, std::vector<double> &a_rhs, double _t0, double _t) {
        double x = a_y[1], y = a_y[2], z = a_y[3];
        double absr = pow(x * x + y * y + z * z, 0.5);
        double dx = (-MU * x) / pow(absr, 3.0 / 2);
        dx -= J2 * MU * Re * Re * 3 * z * z / 2 * ((-5 * x) / pow(absr, 7.0 / 2));
        dx -= MU * Re * Re / 2 * ((-3 * x) / pow(absr, 5.0 / 2));

        double dy = (-MU * y) / pow(absr, 3.0 / 2);
        dy -= J2 * MU * Re * Re * 3 * z * z / 2 * ((-5 * y) / pow(absr, 7.0 / 2));
        dy -= MU * Re * Re / 2 * ((-3 * y) / pow(absr, 5.0 / 2));

        double dz = (-MU * z) / pow(absr, 3.0 / 2);
        dz -= J2 * MU * Re * Re * 3 / 2 * ((2 * x * x * z + 2 * y * y * z - 3 * z * z * z) / pow(absr, 7.0 / 2));
        dz -= MU * Re * Re / 2 * ((-3 * z) / pow(absr, 5.0 / 2));

        a_rhs[0] = dx, a_rhs[1] = dy, a_rhs[2] = dz;

        // and solve from a_y, to a_rhs. Change values of a_rhs
        // a_y -> R ^ {2N + 1}
        // a_rhs -> R ^ {N}

    }
};

template<typename Stepper, typename Observer>
class Everhart_Integrator {
private:
    Stepper *const m_stepper;
    Observer *const m_observer;
public:
    void operator()(double a_t0, const std::vector<double> &a_y0, double a_tEnd, std::vector<double> &a_yEnd) {
        double curr_t = a_t0;
        std::vector<double> curr_y = std::move(a_y0);
        for (int i = 1; curr_t < a_tEnd; i++) {
            std::pair<double, double> res = (*m_stepper)(curr_t, curr_y, 1, a_yEnd);
            curr_t = res.first;
            for (int j = 0; j < a_yEnd.size(); j++) {
                curr_y[j] = a_yEnd[j];
            }
            if (!(*m_observer)(curr_t, i, curr_y)) {
                break;
            }
        }
    }

    Everhart_Integrator(Stepper *a_stepper, Observer *a_observer) : m_stepper(a_stepper), m_observer(a_observer) {}
};

class ObserverPoints {
private:
    std::ofstream ofs;

public:
    explicit ObserverPoints(const std::string &filePath = "../everhart_method/points.txt") : ofs(filePath.c_str(),
                                                                                                 std::ios_base::out) {
    }

    bool operator()(double a_curr_t, int a_n, std::vector<double> &a_curr_y) {
        if (a_curr_t >= 31000000) {
            ofs.close();
            return false; // one year
        }
        ofs <<  a_curr_y[1] << ", " << a_curr_y[2] << ", " << a_curr_y[3] << ",\n"; // x, y, z,
        return true;
    }
};

int main() {
    RHS_U2 u2;
    TimeStepper_everhart stepper(&u2);
    ObserverPoints obs;
    Everhart_Integrator integrator(&stepper, &obs);
    double a = 7500;
    // a = 7500 km (point start is (0, 0, 7500))
    // speed start is (sqrt(mu/a),
    std::vector<double> end(7), start{0, 0, 0, a, pow(MU / a, 0.5), 0, 0};
    integrator(0, start, 31000, end);
    std::cout << end[0] << "\n";
}
