#ifndef FULCRUM_C_H
#define FULCRUM_C_H

#ifdef __cplusplus
extern "C" {
#endif

// 3D Euclidean Vector and Matrix Kernels
void fulcrum_vec3_cross(const double* a, const double* b, double* out);
double fulcrum_vec3_dot(const double* a, const double* b);
void fulcrum_mat3_mul(const double* A, const double* B, double* out);
void fulcrum_mat3_mul_vec(const double* A, const double* x, double* out);

// Spatial 6D Featherstone Motion and Force Transforms
void fulcrum_spatial_transform(const double* X, const double* v, double* out);
void fulcrum_spatial_transform_force(const double* X, const double* f, double* out);
double fulcrum_spatial_kinetic_energy(const double* I, const double* v);

// Contact Mechanics: Projected Gauss-Seidel (PGS) Box-Constrained Solver
// Solves: A * x + b = w, where lo <= x <= hi and complementarity conditions hold
int fulcrum_pgs_solve(
    const double* A,
    const double* b,
    const double* lo,
    const double* hi,
    double* x,
    int n,
    int max_iter,
    double tol
);

// Dense Linear Algebra Accelerators
void fulcrum_dense_mat_mul(
    const double* A,
    const double* B,
    double* C,
    int m,
    int k,
    int n
);

int fulcrum_dense_solve_lu(
    const double* A,
    const double* b,
    double* x,
    int n
);

#ifdef __cplusplus
}
#endif

#endif // FULCRUM_C_H
