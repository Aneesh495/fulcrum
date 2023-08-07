local matrix = require("fulcrum.math.matrix")

local lemke = {}

lemke.STATUS_SOLVED          = "solved"
lemke.STATUS_RAY_TERMINATION = "ray_termination"
lemke.STATUS_MAX_ITER        = "max_iterations"

function lemke.solve(M, q, max_iter, tol)
  -- Solves Linear Complementarity Problem (LCP):
  -- w = M * z + q, subject to w >= 0, z >= 0, w^T * z = 0
  tol = tol or 1e-9
  local n = #q
  max_iter = max_iter or (100 * n)

  -- Check if q >= 0 already (trivial zero solution)
  local all_non_negative = true
  for i = 1, n do
    if q[i] < -tol then
      all_non_negative = false
      break
    end
  end

  if all_non_negative then
    local z_zero = {}
    local w_q = {}
    for i = 1, n do
      z_zero[i] = 0.0
      w_q[i] = math.max(0.0, q[i])
    end
    return lemke.STATUS_SOLVED, z_zero, w_q, 0
  end

  -- Construct initial tableau: [ I, -d, -M | q ]
  -- Variables: w_1..w_n (1..n), z_0 (n+1), z_1..z_n (n+2 .. 2n+1)
  -- Basis initially contains w_1..w_n
  local tab = matrix.zeros(n, 2 * n + 2)

  for i = 1, n do
    tab[i][i] = 1.0           -- w_i column
    tab[i][n + 1] = -1.0      -- -d covering vector (all ones)
    for j = 1, n do
      tab[i][n + 1 + j] = -M[i][j] -- -M columns
    end
    tab[i][2 * n + 2] = q[i]  -- RHS column
  end

  -- Initial basis: basis[i] is index of variable basic in row i
  local basis = {}
  for i = 1, n do
    basis[i] = i -- w_i is in basis
  end

  -- Find row with most negative q_i to pivot z_0 into basis
  local min_q = 0.0
  local pivot_row = 1
  for i = 1, n do
    if tab[i][2 * n + 2] < min_q then
      min_q = tab[i][2 * n + 2]
      pivot_row = i
    end
  end

  -- Pivot z_0 into basis at pivot_row
  local entering_var = n + 1
  local function pivot(row, col)
    local pivot_val = tab[row][col]
    if math.abs(pivot_val) < 1e-15 then
      error("Zero pivot in Lemke algorithm")
    end
    local inv_piv = 1.0 / pivot_val
    for j = 1, 2 * n + 2 do
      tab[row][j] = tab[row][j] * inv_piv
    end

    for i = 1, n do
      if i ~= row then
        local factor = tab[i][col]
        if math.abs(factor) > 1e-15 then
          for j = 1, 2 * n + 2 do
            tab[i][j] = tab[i][j] - factor * tab[row][j]
          end
        end
      end
    end

    basis[row] = col
  end

  local leaving_var = basis[pivot_row]
  pivot(pivot_row, entering_var)

  -- Main Lemke complementary pivoting loop
  local iter = 0
  while iter < max_iter do
    iter = iter + 1

    -- Determine next entering variable (complement of leaving variable)
    if leaving_var <= n then
      -- If w_k left, z_k enters (index n + 1 + k)
      entering_var = n + 1 + leaving_var
    elseif leaving_var > n + 1 then
      -- If z_k left, w_k enters (index k)
      entering_var = leaving_var - (n + 1)
    else
      -- z_0 left: we are done!
      break
    end

    -- Ratio test to determine leaving variable
    local best_ratio = 1e30
    local best_row = nil

    for i = 1, n do
      local col_val = tab[i][entering_var]
      if col_val > tol then
        local rhs_val = tab[i][2 * n + 2]
        local ratio = rhs_val / col_val
        if ratio < best_ratio - 1e-12 then
          best_ratio = ratio
          best_row = i
        elseif math.abs(ratio - best_ratio) <= 1e-12 then
          -- Tie-breaking: Bland rule (choose lower variable index)
          if best_row and basis[i] < basis[best_row] then
            best_row = i
          end
        end
      end
    end

    if not best_row then
      -- Ray termination: no positive entry in pivot column
      return lemke.STATUS_RAY_TERMINATION, nil, nil, iter
    end

    leaving_var = basis[best_row]
    pivot(best_row, entering_var)

    if leaving_var == n + 1 then
      -- z_0 has left basis: optimal complementary solution found
      break
    end
  end

  if iter >= max_iter and leaving_var ~= n + 1 then
    return lemke.STATUS_MAX_ITER, nil, nil, iter
  end

  -- Extract solution: z and w
  local z_sol = {}
  local w_sol = {}
  for i = 1, n do
    z_sol[i] = 0.0
    w_sol[i] = 0.0
  end

  for i = 1, n do
    local var_idx = basis[i]
    local val = math.max(0.0, tab[i][2 * n + 2])
    if var_idx <= n then
      w_sol[var_idx] = val
    elseif var_idx > n + 1 then
      z_sol[var_idx - (n + 1)] = val
    end
  end

  return lemke.STATUS_SOLVED, z_sol, w_sol, iter
end

return lemke
