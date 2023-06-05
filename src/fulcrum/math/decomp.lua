local matrix = require("fulcrum.math.matrix")

local decomp = {}

function decomp.lu(A)
  if A.rows ~= A.cols then
    error("LU decomposition requires square matrix")
  end
  local n = A.rows
  local L = matrix.eye(n)
  local U = matrix.clone(A)
  local P = matrix.eye(n)

  for k = 1, n do
    -- Pivot search
    local max_val = 0.0
    local p = k
    for i = k, n do
      local val = math.abs(U[i][k])
      if val > max_val then
        max_val = val
        p = i
      end
    end

    if max_val < 1e-15 then
      error("Matrix is singular in LU factorisation")
    end

    -- Swap rows in U, P, and preceding columns of L
    if p ~= k then
      U[k], U[p] = U[p], U[k]
      P[k], P[p] = P[p], P[k]
      for j = 1, k - 1 do
        local tmp = L[k][j]
        L[k][j] = L[p][j]
        L[p][j] = tmp
      end
    end

    local piv_inv = 1.0 / U[k][k]
    for i = k + 1, n do
      local factor = U[i][k] * piv_inv
      L[i][k] = factor
      U[i][k] = 0.0
      for j = k + 1, n do
        U[i][j] = U[i][j] - factor * U[k][j]
      end
    end
  end

  return L, U, P
end

function decomp.cholesky(A, regularizer)
  if A.rows ~= A.cols then
    error("Cholesky decomposition requires symmetric square matrix")
  end
  local n = A.rows
  local L = matrix.zeros(n, n)
  local reg = regularizer or 0.0

  for i = 1, n do
    for j = 1, i do
      local sum = 0.0
      for k = 1, j - 1 do
        sum = sum + L[i][k] * L[j][k]
      end

      if i == j then
        local diag = A[i][i] + reg - sum
        if diag <= 1e-15 then
          error(string.format("Matrix is not positive-definite at index %d (diag = %g)", i, diag))
        end
        L[i][j] = math.sqrt(diag)
      else
        L[i][j] = (A[i][j] - sum) / L[j][j]
      end
    end
  end

  return L
end

function decomp.cholesky_solve(L, b)
  local n = L.rows
  -- Solve L y = b (forward substitution)
  local y = {}
  for i = 1, n do
    local sum = b[i]
    for j = 1, i - 1 do
      sum = sum - L[i][j] * y[j]
    end
    y[i] = sum / L[i][i]
  end

  -- Solve L^T x = y (backward substitution)
  local x = {}
  for i = n, 1, -1 do
    local sum = y[i]
    for j = i + 1, n do
      sum = sum - L[j][i] * x[j]
    end
    x[i] = sum / L[i][i]
  end

  return x
end

function decomp.qr(A)
  local m = A.rows
  local n = A.cols
  local Q = matrix.eye(m)
  local R = matrix.clone(A)

  local num_reflections = math.min(m - 1, n)
  for k = 1, num_reflections do
    -- Extract column vector below diagonal
    local norm_x_sq = 0.0
    for i = k, m do
      local v = R[i][k]
      norm_x_sq = norm_x_sq + v * v
    end
    local norm_x = math.sqrt(norm_x_sq)

    if norm_x > 1e-15 then
      local alpha = (R[k][k] >= 0) and (-norm_x) or norm_x
      local u1 = R[k][k] - alpha
      local norm_u = math.sqrt(norm_x_sq - R[k][k] * R[k][k] + u1 * u1)
      local inv_norm_u = 1.0 / norm_u

      local v = {}
      for i = 1, m do
        if i < k then
          v[i] = 0.0
        elseif i == k then
          v[i] = u1 * inv_norm_u
        else
          v[i] = R[i][k] * inv_norm_u
        end
      end

      -- Apply Householder transformation to R: R = (I - 2 v v^T) R
      for j = k, n do
        local dot = 0.0
        for i = k, m do
          dot = dot + v[i] * R[i][j]
        end
        local tau = 2.0 * dot
        for i = k, m do
          R[i][j] = R[i][j] - tau * v[i]
        end
      end

      -- Accumulate Q: Q = Q (I - 2 v v^T)
      for i = 1, m do
        local dot = 0.0
        for p = k, m do
          dot = dot + Q[i][p] * v[p]
        end
        local tau = 2.0 * dot
        for p = k, m do
          Q[i][p] = Q[i][p] - tau * v[p]
        end
      end
    end
  end

  return Q, R
end

function decomp.qr_solve(Q, R, b)
  -- Solve R x = Q^T b
  local m = R.rows
  local n = R.cols
  local Qt = matrix.transpose(Q)
  local Qt_b = matrix.mul_vec(Qt, b)

  local x = {}
  for i = n, 1, -1 do
    local sum = Qt_b[i]
    for j = i + 1, n do
      sum = sum - R[i][j] * x[j]
    end
    if math.abs(R[i][i]) < 1e-15 then
      x[i] = 0.0
    else
      x[i] = sum / R[i][i]
    end
  end

  return x
end

function decomp.svd(A, max_iter, tol)
  -- Two-sided Jacobi rotation SVD for m >= n
  local m = A.rows
  local n = A.cols
  local swapped = false
  local B = A

  if m < n then
    B = matrix.transpose(A)
    m, n = n, m
    swapped = true
  end

  max_iter = max_iter or 100
  tol = tol or 1e-12

  local U = matrix.clone(B)
  local V = matrix.eye(n)
  local S = {}

  for _ = 1, max_iter do
    local max_off = 0.0
    for i = 1, n - 1 do
      for j = i + 1, n do
        -- Compute Gram matrix elements (U^T U)
        local a = 0.0
        local b_elem = 0.0
        local c = 0.0
        for k = 1, m do
          local uki = U[k][i]
          local ukj = U[k][j]
          a = a + uki * uki
          b_elem = b_elem + ukj * ukj
          c = c + uki * ukj
        end

        max_off = math.max(max_off, math.abs(c))

        if math.abs(c) > tol then
          local zeta = (b_elem - a) / (2.0 * c)
          local t
          if zeta >= 0 then
            t = 1.0 / (zeta + math.sqrt(1.0 + zeta * zeta))
          else
            t = -1.0 / (-zeta + math.sqrt(1.0 + zeta * zeta))
          end
          local cos_theta = 1.0 / math.sqrt(1.0 + t * t)
          local sin_theta = t * cos_theta

          -- Rotate columns of U
          for k = 1, m do
            local uki = U[k][i]
            local ukj = U[k][j]
            U[k][i] = cos_theta * uki - sin_theta * ukj
            U[k][j] = sin_theta * uki + cos_theta * ukj
          end

          -- Rotate columns of V
          for k = 1, n do
            local vki = V[k][i]
            local vkj = V[k][j]
            V[k][i] = cos_theta * vki - sin_theta * vkj
            V[k][j] = sin_theta * vki + cos_theta * vkj
          end
        end
      end
    end

    if max_off < tol then break end
  end

  -- Compute singular values and normalise U
  for j = 1, n do
    local col_norm_sq = 0.0
    for k = 1, m do
      col_norm_sq = col_norm_sq + U[k][j] * U[k][j]
    end
    local s = math.sqrt(col_norm_sq)
    S[j] = s
    if s > 1e-15 then
      local inv_s = 1.0 / s
      for k = 1, m do
        U[k][j] = U[k][j] * inv_s
      end
    end
  end

  if swapped then
    -- A = (U S V^T)^T = V S U^T
    return V, S, U
  else
    return U, S, V
  end
end

function decomp.pinv(A, tol)
  local U, S, V = decomp.svd(A)
  tol = tol or 1e-10

  -- A^+ = V * S^+ * U^T
  local m = A.rows
  local n = A.cols
  local k = #S
  local pinv_mat = matrix.zeros(n, m)

  for i = 1, n do
    for j = 1, m do
      local sum = 0.0
      for p = 1, k do
        if S[p] > tol then
          sum = sum + V[i][p] * (1.0 / S[p]) * U[j][p]
        end
      end
      pinv_mat[i][j] = sum
    end
  end

  return pinv_mat
end

return decomp
