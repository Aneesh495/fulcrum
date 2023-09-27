local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")

local c_accel = {}

local has_ffi, ffi = pcall(require, "ffi")
local clib = nil

if has_ffi then
  ffi.cdef[[
    void fulcrum_vec3_cross(const double* a, const double* b, double* out);
    double fulcrum_vec3_dot(const double* a, const double* b);
    void fulcrum_mat3_mul(const double* A, const double* B, double* out);
    void fulcrum_mat3_mul_vec(const double* A, const double* x, double* out);
    void fulcrum_spatial_transform(const double* X, const double* v, double* out);
    void fulcrum_spatial_transform_force(const double* X, const double* f, double* out);
    double fulcrum_spatial_kinetic_energy(const double* I, const double* v);
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
  ]]

  local candidate_paths = {
    "./src/c/libfulcrum.dylib",
    "./src/c/libfulcrum.so",
    "src/c/libfulcrum.dylib",
    "src/c/libfulcrum.so"
  }

  for _, path in ipairs(candidate_paths) do
    local ok, lib = pcall(ffi.load, path)
    if ok and lib then
      clib = lib
      break
    end
  end
end

function c_accel.is_available()
  return clib ~= nil
end

function c_accel.vec3_cross(a, b)
  if clib then
    local a_arr = ffi.new("double[3]", { a.x, a.y, a.z })
    local b_arr = ffi.new("double[3]", { b.x, b.y, b.z })
    local out = ffi.new("double[3]")
    clib.fulcrum_vec3_cross(a_arr, b_arr, out)
    return vec3.new(out[0], out[1], out[2])
  else
    return vec3.cross(a, b)
  end
end

function c_accel.vec3_dot(a, b)
  if clib then
    local a_arr = ffi.new("double[3]", { a.x, a.y, a.z })
    local b_arr = ffi.new("double[3]", { b.x, b.y, b.z })
    return clib.fulcrum_vec3_dot(a_arr, b_arr)
  else
    return vec3.dot(a, b)
  end
end

function c_accel.mat3_mul(A, B)
  if clib then
    local A_arr = ffi.new("double[9]", {
      A[1][1], A[1][2], A[1][3],
      A[2][1], A[2][2], A[2][3],
      A[3][1], A[3][2], A[3][3]
    })
    local B_arr = ffi.new("double[9]", {
      B[1][1], B[1][2], B[1][3],
      B[2][1], B[2][2], B[2][3],
      B[3][1], B[3][2], B[3][3]
    })
    local out = ffi.new("double[9]")
    clib.fulcrum_mat3_mul(A_arr, B_arr, out)
    return mat3.new(
      out[0], out[1], out[2],
      out[3], out[4], out[5],
      out[6], out[7], out[8]
    )
  else
    return mat3.mul(A, B)
  end
end

function c_accel.spatial_transform(X, v)
  if clib then
    local X_arr = ffi.new("double[36]")
    for i = 1, 6 do
      for j = 1, 6 do
        X_arr[(i - 1) * 6 + (j - 1)] = X[i][j]
      end
    end
    local v_arr = ffi.new("double[6]", v)
    local out = ffi.new("double[6]")
    clib.fulcrum_spatial_transform(X_arr, v_arr, out)
    local res = {}
    for i = 1, 6 do res[i] = out[i - 1] end
    return res
  else
    return matrix.mul_vec(X, v)
  end
end

function c_accel.spatial_transform_force(X, f)
  if clib then
    local X_arr = ffi.new("double[36]")
    for i = 1, 6 do
      for j = 1, 6 do
        X_arr[(i - 1) * 6 + (j - 1)] = X[i][j]
      end
    end
    local f_arr = ffi.new("double[6]", f)
    local out = ffi.new("double[6]")
    clib.fulcrum_spatial_transform_force(X_arr, f_arr, out)
    local res = {}
    for i = 1, 6 do res[i] = out[i - 1] end
    return res
  else
    local Xt = matrix.transpose(X)
    return matrix.mul_vec(Xt, f)
  end
end

function c_accel.spatial_kinetic_energy(I_sp, v)
  if clib then
    local I_arr = ffi.new("double[36]")
    for i = 1, 6 do
      for j = 1, 6 do
        I_arr[(i - 1) * 6 + (j - 1)] = I_sp[i][j]
      end
    end
    local v_arr = ffi.new("double[6]", v)
    return clib.fulcrum_spatial_kinetic_energy(I_arr, v_arr)
  else
    local temp = matrix.mul_vec(I_sp, v)
    local sum = 0.0
    for i = 1, 6 do sum = sum + v[i] * temp[i] end
    return 0.5 * sum
  end
end

function c_accel.pgs_solve(A, b, lo, hi, x_init, max_iter, tol)
  local n = #b
  max_iter = max_iter or 100
  tol = tol or 1e-6

  if clib then
    local A_arr = ffi.new("double[?]", n * n)
    for i = 1, n do
      for j = 1, n do
        A_arr[(i - 1) * n + (j - 1)] = A[i][j]
      end
    end
    local b_arr = ffi.new("double[?]", n, b)
    local x_arr = ffi.new("double[?]", n)
    for i = 1, n do
      x_arr[i - 1] = x_init and x_init[i] or 0.0
    end

    local lo_arr = lo and ffi.new("double[?]", n, lo) or nil
    local hi_arr = hi and ffi.new("double[?]", n, hi) or nil

    local iters = clib.fulcrum_pgs_solve(A_arr, b_arr, lo_arr, hi_arr, x_arr, n, max_iter, tol)

    local res = {}
    for i = 1, n do res[i] = x_arr[i - 1] end
    return res, iters
  else
    -- Pure Lua Projected Gauss-Seidel fallback
    local x = {}
    for i = 1, n do x[i] = x_init and x_init[i] or 0.0 end

    for iter = 1, max_iter do
      local max_diff = 0.0
      for i = 1, n do
        local sum = b[i]
        local row = A[i]
        for j = 1, n do
          if j ~= i then
            sum = sum + row[j] * x[j]
          end
        end

        local diag = row[i]
        if math.abs(diag) < 1e-15 then
          diag = (diag >= 0) and 1e-15 or -1e-15
        end

        local new_val = -sum / diag
        if lo and new_val < lo[i] then new_val = lo[i] end
        if hi and new_val > hi[i] then new_val = hi[i] end

        local diff = math.abs(new_val - x[i])
        if diff > max_diff then max_diff = diff end
        x[i] = new_val
      end

      if max_diff < tol then
        return x, iter
      end
    end

    return x, max_iter
  end
end

return c_accel
