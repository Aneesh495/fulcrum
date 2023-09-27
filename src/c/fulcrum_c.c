#include "fulcrum_c.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void fulcrum_vec3_cross(const double* a, const double* b, double* out) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

double fulcrum_vec3_dot(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void fulcrum_mat3_mul(const double* A, const double* B, double* out) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += A[i * 3 + k] * B[k * 3 + j];
            }
            out[i * 3 + j] = sum;
        }
    }
}

void fulcrum_mat3_mul_vec(const double* A, const double* x, double* out) {
    for (int i = 0; i < 3; ++i) {
        out[i] = A[i * 3 + 0] * x[0] +
                 A[i * 3 + 1] * x[1] +
                 A[i * 3 + 2] * x[2];
    }
}

void fulcrum_spatial_transform(const double* X, const double* v, double* out) {
    for (int i = 0; i < 6; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 6; ++j) {
            sum += X[i * 6 + j] * v[j];
        }
        out[i] = sum;
    }
}

void fulcrum_spatial_transform_force(const double* X, const double* f, double* out) {
    // Force transforms via transpose of motion transform matrix (power invariance)
    for (int i = 0; i < 6; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 6; ++j) {
            sum += X[j * 6 + i] * f[j];
        }
        out[i] = sum;
    }
}

double fulcrum_spatial_kinetic_energy(const double* I, const double* v) {
    double temp[6];
    for (int i = 0; i < 6; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 6; ++j) {
            sum += I[i * 6 + j] * v[j];
        }
        temp[i] = sum;
    }

    double energy = 0.0;
    for (int i = 0; i < 6; ++i) {
        energy += v[i] * temp[i];
    }
    return 0.5 * energy;
}

int fulcrum_pgs_solve(
    const double* A,
    const double* b,
    const double* lo,
    const double* hi,
    double* x,
    int n,
    int max_iter,
    double tol
) {
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            const double* row = &A[i * n];

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += row[j] * x[j];
                }
            }

            double diag = row[i];
            if (fabs(diag) < 1e-15) {
                diag = (diag >= 0.0) ? 1e-15 : -1e-15;
            }

            double new_val = -sum / diag;
            if (lo && new_val < lo[i]) new_val = lo[i];
            if (hi && new_val > hi[i]) new_val = hi[i];

            double diff = fabs(new_val - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }

            x[i] = new_val;
        }

        if (max_diff < tol) {
            return iter + 1; // Converged
        }
    }

    return max_iter;
}

void fulcrum_dense_mat_mul(
    const double* A,
    const double* B,
    double* C,
    int m,
    int k,
    int n
) {
    memset(C, 0, sizeof(double) * m * n);
    for (int i = 0; i < m; ++i) {
        for (int p = 0; p < k; ++p) {
            double a_ip = A[i * k + p];
            if (a_ip != 0.0) {
                for (int j = 0; j < n; ++j) {
                    C[i * n + j] += a_ip * B[p * n + j];
                }
            }
        }
    }
}

int fulcrum_dense_solve_lu(
    const double* A,
    const double* b,
    double* x,
    int n
) {
    double* LU = (double*)malloc(sizeof(double) * n * n);
    int* p = (int*)malloc(sizeof(int) * n);
    if (!LU || !p) {
        free(LU);
        free(p);
        return 0;
    }
    memcpy(LU, A, sizeof(double) * n * n);

    for (int i = 0; i < n; ++i) p[i] = i;

    // LU factorisation with partial pivoting
    for (int k = 0; k < n; ++k) {
        double max_val = 0.0;
        int pivot = k;
        for (int i = k; i < n; ++i) {
            double val = fabs(LU[i * n + k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        if (max_val < 1e-15) {
            free(LU);
            free(p);
            return 0; // Singular matrix
        }

        if (pivot != k) {
            int tmp_p = p[k];
            p[k] = p[pivot];
            p[pivot] = tmp_p;

            for (int j = 0; j < n; ++j) {
                double tmp = LU[k * n + j];
                LU[k * n + j] = LU[pivot * n + j];
                LU[pivot * n + j] = tmp;
            }
        }

        double inv_diag = 1.0 / LU[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            LU[i * n + k] *= inv_diag;
            double factor = LU[i * n + k];
            for (int j = k + 1; j < n; ++j) {
                LU[i * n + j] -= factor * LU[k * n + j];
            }
        }
    }

    // Forward solve: L y = P b
    double* y = (double*)malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i) {
        double sum = b[p[i]];
        for (int j = 0; j < i; ++j) {
            sum -= LU[i * n + j] * y[j];
        }
        y[i] = sum;
    }

    // Backward solve: U x = y
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= LU[i * n + j] * x[j];
        }
        x[i] = sum / LU[i * n + i];
    }

    free(y);
    free(LU);
    free(p);
    return 1; // Success
}
