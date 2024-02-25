#include <gsl/gsl_linalg.h>

void solve(int N) {
    auto *a_data = (double *) calloc(N * N, sizeof(double));
    auto *b_data = (double *) calloc(N, sizeof(double));
    if (!a_data || !b_data) {
        return;
    }

    for (int k = 0; k < N - 1; ++k) {
        a_data[N * k + k] = -1;
        for (int n = k + 1; n < N; ++n) {
            if (n % 2 == 0) {
                if (k == 1) {
                    a_data[N * k + n] = n;
                } else if (k % 2 == 1) {
                    a_data[N * k + n] = 2 * n;
                }
            } else {
                if (k == 0) {
                    a_data[N * k + n] = n;
                } else if (k % 2 == 0) {
                    a_data[N * k + n] = 2 * n;
                }
            }
        }
    }
    for (int i = 0; i < N; i++) {
        double val = 0;
        if (i % 4 == 0) val = 1;
        if (i % 4 == 2) val = -1;
        a_data[N * (N - 1) + i] = val;
    }
    b_data[N - 1] = 1;

    gsl_matrix_view m = gsl_matrix_view_array(a_data, N, N);
    gsl_vector_view b = gsl_vector_view_array(b_data, N);
    gsl_vector *x = gsl_vector_alloc(N);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

    gsl_vector_fprintf(stdout, x, "%.17g");

    gsl_permutation_free(p);
    gsl_vector_free(x);
    free(b_data);
    free(a_data);
}

int main() {
    solve(15);
    return 0;
}
