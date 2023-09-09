local matrix = require("fulcrum.math.matrix")

local collocation = {}

function collocation.hermite_simpson_defect(x_k, u_k, x_next, u_next, dynamics, h)
  local f_k = dynamics(x_k, u_k)
  local f_next = dynamics(x_next, u_next)

  local nx = #x_k
  local nu = #u_k

  -- Midpoint state: x_mid = 0.5*(x_k + x_next) + (h/8)*(f_k - f_next)
  local x_mid = {}
  for i = 1, nx do
    x_mid[i] = 0.5 * (x_k[i] + x_next[i]) + (h / 8.0) * (f_k[i] - f_next[i])
  end

  -- Midpoint control: u_mid = 0.5*(u_k + u_next)
  local u_mid = {}
  for i = 1, nu do
    u_mid[i] = 0.5 * (u_k[i] + u_next[i])
  end

  local f_mid = dynamics(x_mid, u_mid)

  -- Defect constraint: x_next - x_k - (h/6)*(f_k + 4*f_mid + f_next) = 0
  local defect = {}
  for i = 1, nx do
    defect[i] = x_next[i] - x_k[i] - (h / 6.0) * (f_k[i] + 4.0 * f_mid[i] + f_next[i])
  end

  return defect, x_mid, u_mid, f_k, f_mid, f_next
end

function collocation.solve(prob, opts)
  opts = opts or {}
  local M = prob.segments or 8
  local dt = prob.duration / M
  local nx = prob.n_x
  local nu = prob.n_u

  local max_iter = opts.max_iter or 15
  local tol_defect = opts.tol_defect or 1e-4
  local reg_kkt = opts.reg or 1e-6

  local n_vars = (M + 1) * (nx + nu)

  local function unpack_z(Z)
    local x_traj = {}
    local u_traj = {}
    local idx = 1
    for k = 1, M + 1 do
      local x = {}
      for i = 1, nx do
        x[i] = Z[idx]
        idx = idx + 1
      end
      local u = {}
      for i = 1, nu do
        u[i] = Z[idx]
        idx = idx + 1
      end
      x_traj[k] = x
      u_traj[k] = u
    end
    return x_traj, u_traj
  end

  local function pack_z(x_traj, u_traj)
    local Z = {}
    for k = 1, M + 1 do
      for i = 1, nx do table.insert(Z, x_traj[k][i]) end
      for i = 1, nu do table.insert(Z, u_traj[k][i]) end
    end
    return Z
  end

  local function evaluate_constraints(Z)
    local x_t, u_t = unpack_z(Z)
    local c = {}

    -- Initial state condition: x(0) - x_start = 0
    for i = 1, nx do
      table.insert(c, x_t[1][i] - prob.x_start[i])
    end

    -- Defect constraints across segments
    for k = 1, M do
      local def = collocation.hermite_simpson_defect(
        x_t[k], u_t[k], x_t[k + 1], u_t[k + 1], prob.dynamics, dt
      )
      for i = 1, nx do
        table.insert(c, def[i])
      end
    end

    -- Terminal state condition: x(M+1) - x_goal = 0
    if prob.x_goal then
      for i = 1, nx do
        table.insert(c, x_t[M + 1][i] - prob.x_goal[i])
      end
    end

    return c
  end

  local function evaluate_cost(Z)
    local x_t, u_t = unpack_z(Z)
    local total = 0.0

    for k = 1, M do
      local _, x_mid, u_mid = collocation.hermite_simpson_defect(
        x_t[k], u_t[k], x_t[k + 1], u_t[k + 1], prob.dynamics, dt
      )
      local ck = prob.cost and prob.cost(x_t[k], u_t[k]) or 0.0
      local cm = prob.cost and prob.cost(x_mid, u_mid) or 0.0
      local cn = prob.cost and prob.cost(x_t[k + 1], u_t[k + 1]) or 0.0
      total = total + (dt / 6.0) * (ck + 4.0 * cm + cn)
    end

    if prob.terminal_cost then
      total = total + prob.terminal_cost(x_t[M + 1])
    end

    return total
  end

  -- Initial linear interpolation trajectory
  local x_init = {}
  local u_init = {}
  for k = 1, M + 1 do
    local alpha = (k - 1) / M
    local x_k = {}
    for i = 1, nx do
      local s = prob.x_start[i] or 0.0
      local g = prob.x_goal and prob.x_goal[i] or s
      x_k[i] = (1.0 - alpha) * s + alpha * g
    end
    x_init[k] = x_k

    local u_k = {}
    for i = 1, nu do u_k[i] = 0.0 end
    u_init[k] = u_k
  end

  local Z = pack_z(x_init, u_init)
  local converged = false
  local final_defect = 0.0
  local final_cost = 0.0

  for _ = 1, max_iter do
    local c = evaluate_constraints(Z)
    local max_c = 0.0
    for _, v in ipairs(c) do
      local abs_v = math.abs(v)
      if abs_v > max_c then max_c = abs_v end
    end

    final_defect = max_c
    final_cost = evaluate_cost(Z)

    if max_c < tol_defect then
      converged = true
      break
    end

    -- Numerical differentiation for constraint Jacobian and cost gradient
    local eps = 1e-6
    local inv_eps = 1.0 / eps
    local n_c = #c
    local Jc = matrix.zeros(n_c, n_vars)
    local g = {}

    for j = 1, n_vars do
      local orig = Z[j]
      Z[j] = orig + eps
      local c_plus = evaluate_constraints(Z)
      local cost_plus = evaluate_cost(Z)
      Z[j] = orig

      for i = 1, n_c do
        Jc[i][j] = (c_plus[i] - c[i]) * inv_eps
      end
      g[j] = (cost_plus - final_cost) * inv_eps
    end

    -- Form KKT system: [ H, Jc^T ; Jc, -reg*I ] * [ dZ ; lambda ] = [ -g ; -c ]
    local kkt_dim = n_vars + n_c
    local KKT = matrix.zeros(kkt_dim, kkt_dim)
    local rhs = {}

    -- Cost Hessian block (regularized for positive definiteness)
    for i = 1, n_vars do
      KKT[i][i] = 1.0
      rhs[i] = -g[i]
    end

    for i = 1, n_c do
      rhs[n_vars + i] = -c[i]
      for j = 1, n_vars do
        local jac_val = Jc[i][j]
        KKT[n_vars + i][j] = jac_val
        KKT[j][n_vars + i] = jac_val
      end
      KKT[n_vars + i][n_vars + i] = -reg_kkt
    end

    local ok, sol = pcall(function() return matrix.solve_lu(KKT, rhs) end)
    if not ok then
      break
    end

    -- Backtracking line search on merit function: cost + 100 * ||c||
    local alpha = 1.0
    local base_merit = final_cost + 100.0 * max_c
    local accepted = false

    for _ = 1, 6 do
      local cand_Z = {}
      for j = 1, n_vars do
        cand_Z[j] = Z[j] + alpha * sol[j]
      end

      local cand_c = evaluate_constraints(cand_Z)
      local cand_max_c = 0.0
      for _, v in ipairs(cand_c) do
        local abs_v = math.abs(v)
        if abs_v > cand_max_c then cand_max_c = abs_v end
      end
      local cand_cost = evaluate_cost(cand_Z)
      local cand_merit = cand_cost + 100.0 * cand_max_c

      if cand_merit < base_merit or cand_max_c < max_c then
        Z = cand_Z
        accepted = true
        break
      end
      alpha = alpha * 0.5
    end

    if not accepted then
      -- Fallback step with small relaxation
      for j = 1, n_vars do
        Z[j] = Z[j] + 0.1 * sol[j]
      end
    end
  end

  local x_res, u_res = unpack_z(Z)

  return {
    x = x_res,
    u = u_res,
    M = M,
    dt = dt,
    cost = final_cost,
    defect_norm = final_defect,
    converged = converged
  }
end

return collocation
