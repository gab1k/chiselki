#include <iostream>
#include <vector>

std::vector<int> get_c(int N) {
    std::vector<std::vector<int>> a(N + 1, std::vector<int>(N + 1, 0));
    std::vector<int> T(N + 1, 0);
    for (int n = 0; n < N + 1; n++) {
        for (int k = 0; k < N + 1; k++) {
            if (n & 1) { // n = odd
                a[n][k] = (k & 1) * (k != 0 ? 2 * n : n);
            } else {
                a[n][k] = (1 - k & 1) * (k != 1 ? 2 * n : n);
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