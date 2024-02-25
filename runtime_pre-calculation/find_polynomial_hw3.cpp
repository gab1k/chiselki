#include <iostream>
#include <vector>

std::vector<int> get_c(int N) {
    std::vector<std::vector<int>> a(N + 1, std::vector<int>(N + 1, 0));
    std::vector<int> T(N + 1, 0);
    for (int n = 0; n < N + 1; n++) {
        for (int k = 0; k < N + 1; k++) {
            if (n & 1) { // n = odd
                if (k & 1) { // k = odd
                    a[n][k] = 0;
                } else if (k == 0) {
                    a[n][k] = n;
                } else {
                    a[n][k] = 2 * n;
                }
            } else {
                if (k & 1) {
                    a[n][k] = 2 * n;
                }
            }
        }
    }
    for (int n = 0; n < N + 1; n += 2) {
        if ((n / 2) & 1) { // n/2 odd
            T[n] = -1;
        } else {
            T[n] = 1;
        }
    }
    return {};
}

int main() {
    std::vector<int> c = get_c(17);
    for (int &c_i: c) {
        std::cout << c_i << " ";
    }
}