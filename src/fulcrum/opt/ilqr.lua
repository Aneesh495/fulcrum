local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")

local ilqr = {}

local function finite_diff_derivatives(prob, x, u, k)
  local nx = prob.n_x
  local nu = prob.n_u
  local h = 1e-6
  local inv_2h = 1.0 / (2.0 * h)

  -- Dynamics Jacobians A = df/dx, B = df/du
  local A = matrix.zeros(nx, nx)
  local B = matrix.zeros(nx, nu)

  for i = 1, nx do
    local x_plus = {}
    local x_minus = {}
    for j = 1, nx do
      x_plus[j] = x[j] + (i == j and h or 0.0)
      x_minus[j] = x[j] - (i == j and h or 0.0)
    end
    local f_plus = prob.dynamics(x_plus, u, k)
    local f_minus = prob.dynamics(x_minus, u, k)
    for r = 1, nx do
      A[r][i] = (f_plus[r] - f_minus[r]) * inv_2h
    end
  end

  for i = 1, nu do
    local u_plus = {}
    local u_minus = {}
    for j = 1, nu do
      u_plus[j] = u[j] + (i == j and h or 0.0)
      u_minus[j] = u[j] - (i == j and h or 0.0)
    end
    local f_plus = prob.dynamics(x, u_plus, k)
    local f_minus = prob.dynamics(x, u_minus, k)
    for r = 1, nx do
      B[r][i] = (f_plus[r] - f_minus[r]) * inv_2h
    end
  end

  -- Cost gradient lx, lu
  local lx = {}
  for i = 1, nx do
    local x_plus = {}
    local x_minus = {}
    for j = 1, nx do
      x_plus[j] = x[j] + (i == j and h or 0.0)
      x_minus[j] = x[j] - (i == j and h or 0.0)
    end
    lx[i] = (prob.cost(x_plus, u, k) - prob.cost(x_minus, u, k)) * inv_2h
  end

  local lu = {}
  for i = 1, nu do
    local u_plus = {}
    local u_minus = {}
    for j = 1, nu do
      u_plus[j] = u[j] + (i == j and h or 0.0)
      u_minus[j] = u[j] - (i == j and h or 0.0)
    end
    lu[i] = (prob.cost(x, u_plus, k) - prob.cost(x, u_minus, k)) * inv_2h
  end

  -- Cost Hessian lxx, luu, lux
  local lxx = matrix.zeros(nx, nx)
  local c_base = prob.cost(x, u, k)
  local inv_h2 = 1.0 / (h * h)

  for i = 1, nx do
    for j = i, nx do
      if i == j then
        local xp, xm = {}, {}
        for m = 1, nx do
          xp[m] = x[m] + (m == i and h or 0.0)
          xm[m] = x[m] - (m == i and h or 0.0)
        end
        local val = (prob.cost(xp, u, k) - 2.0 * c_base + prob.cost(xm, u, k)) * inv_h2
        lxx[i][i] = val
      else
        local xpp, xpm, xmp, xmm = {}, {}, {}, {}
        for m = 1, nx do
          local di = (m == i and h or 0.0)
          local dj = (m == j and h or 0.0)
          xpp[m] = x[m] + di + dj
          xpm[m] = x[m] + di - dj
          xmp[m] = x[m] - di + dj
          xmm[m] = x[m] - di - dj
        end
        local val = (prob.cost(xpp, u, k) - prob.cost(xpm, u, k) - prob.cost(xmp, u, k) + prob.cost(xmm, u, k)) * (0.25 * inv_h2)
        lxx[i][j] = val
        lxx[j][i] = val
      end
    end
  end

  local luu = matrix.zeros(nu, nu)
  for i = 1, nu do
    for j = i, nu do
      if i == j then
        local up, um = {}, {}
        for m = 1, nu do
          up[m] = u[m] + (m == i and h or 0.0)
          um[m] = u[m] - (m == i and h or 0.0)
        end
        local val = (prob.cost(x, up, k) - 2.0 * c_base + prob.cost(x, um, k)) * inv_h2
        luu[i][i] = val
      else
        local upp, upm, ump, umm = {}, {}, {}, {}
        for m = 1, nu do
          local di = (m == i and h or 0.0)
          local dj = (m == j and h or 0.0)
          upp[m] = u[m] + di + dj
          upm[m] = u[m] + di - dj
          ump[m] = u[m] - di + dj
          umm[m] = u[m] - di - dj
        end
        local val = (prob.cost(x, upp, k) - prob.cost(x, upm, k) - prob.cost(x, ump, k) + prob.cost(x, umm, k)) * (0.25 * inv_h2)
        luu[i][j] = val
        luu[j][i] = val
      end
    end
  end

  local lux = matrix.zeros(nu, nx)
  for i = 1, nu do
    for j = 1, nx do
      local xp_up, xp_um, xm_up, xm_um = {}, {}, {}, {}
      local up, um = {}, {}
      for m = 1, nu do
        up[m] = u[m] + (m == i and h or 0.0)
        um[m] = u[m] - (m == i and h or 0.0)
      end
      local xp, xm = {}, {}
      for m = 1, nx do
        xp[m] = x[m] + (m == j and h or 0.0)
        xm[m] = x[m] - (m == j and h or 0.0)
      end
      local val = (prob.cost(xp, up, k) - prob.cost(xp, um, k) - prob.cost(xm, up, k) + prob.cost(xm, um, k)) * (0.25 * inv_h2)
      lux[i][j] = val
    end
  end

  return A, B, lx, lu, lxx, luu, lux
end

local function finite_diff_terminal(prob, x)
  local nx = prob.n_x
  local h = 1e-6
  local inv_2h = 1.0 / (2.0 * h)
  local inv_h2 = 1.0 / (h * h)
  local c_base = prob.terminal_cost(x)

  local Vx = {}
  for i = 1, nx do
    local xp, xm = {}, {}
    for m = 1, nx do
      xp[m] = x[m] + (m == i and h or 0.0)
      xm[m] = x[m] - (m == i and h or 0.0)
    end
    Vx[i] = (prob.terminal_cost(xp) - prob.terminal_cost(xm)) * inv_2h
  end

  local Vxx = matrix.zeros(nx, nx)
  for i = 1, nx do
    for j = i, nx do
      if i == j then
        local xp, xm = {}, {}
        for m = 1, nx do
          xp[m] = x[m] + (m == i and h or 0.0)
          xm[m] = x[m] - (m == i and h or 0.0)
        end
        local val = (prob.terminal_cost(xp) - 2.0 * c_base + prob.terminal_cost(xm)) * inv_h2
        Vxx[i][i] = val
      else
        local xpp, xpm, xmp, xmm = {}, {}, {}, {}
        for m = 1, nx do
          local di = (m == i and h or 0.0)
          local dj = (m == j and h or 0.0)
          xpp[m] = x[m] + di + dj
          xpm[m] = x[m] + di - dj
          xmp[m] = x[m] - di + dj
          xmm[m] = x[m] - di - dj
        end
        local val = (prob.terminal_cost(xpp) - prob.terminal_cost(xpm) - prob.terminal_cost(xmp) + prob.terminal_cost(xmm)) * (0.25 * inv_h2)
        Vxx[i][j] = val
        Vxx[j][i] = val
      end
    end
  end

  return Vx, Vxx
end

function ilqr.rollout(prob, x0, u_traj)
  local N = prob.N
  local x_traj = { x0 }
  local total_cost = 0.0

  for k = 0, N - 1 do
    local u = u_traj[k + 1]
    local x = x_traj[k + 1]
    total_cost = total_cost + prob.cost(x, u, k)
    local x_next = prob.dynamics(x, u, k)
    table.insert(x_traj, x_next)
  end

  total_cost = total_cost + prob.terminal_cost(x_traj[N + 1])
  return x_traj, total_cost
end

function ilqr.solve(prob, x0, u_init, opts)
  opts = opts or {}
  local max_iter = opts.max_iter or 40
  local tol = opts.tol or 1e-4
  local mu = opts.mu_init or 1e-6
  local mu_min = opts.mu_min or 1e-8
  local mu_max = opts.mu_max or 1e8
  local delta_mu = opts.delta_mu or 2.0
  local alphas = opts.alphas or { 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125 }

  local N = prob.N
  local nx = prob.n_x
  local nu = prob.n_u

  -- Initial control trajectory
  local u_traj = {}
  for k = 1, N do
    if u_init and u_init[k] then
      u_traj[k] = u_init[k]
    else
      local u_zero = {}
      for j = 1, nu do u_zero[j] = 0.0 end
      u_traj[k] = u_zero
    end
  end

  local x_traj, current_cost = ilqr.rollout(prob, x0, u_traj)
  local cost_history = { current_cost }
  local K_traj = {}
  local k_traj = {}

  local converged = false
  local iter_count = 0

  for iter = 1, max_iter do
    iter_count = iter
    -- Backward Pass
    local Vx, Vxx
    if prob.terminal_derivatives then
      Vx, Vxx = prob.terminal_derivatives(x_traj[N + 1])
    else
      Vx, Vxx = finite_diff_terminal(prob, x_traj[N + 1])
    end

    local backward_success = true
    local expected_cost_reduction = 0.0

    for k = N - 1, 0, -1 do
      local x = x_traj[k + 1]
      local u = u_traj[k + 1]

      local A, B, lx, lu, lxx, luu, lux
      if prob.derivatives then
        A, B, lx, lu, lxx, luu, lux = prob.derivatives(x, u, k)
      else
        A, B, lx, lu, lxx, luu, lux = finite_diff_derivatives(prob, x, u, k)
      end

      -- Qx = lx + A^T * Vx
      local At = matrix.transpose(A)
      local Bt = matrix.transpose(B)
      local At_Vx = matrix.mul_vec(At, Vx)
      local Bt_Vx = matrix.mul_vec(Bt, Vx)

      local Qx = {}
      for i = 1, nx do Qx[i] = lx[i] + At_Vx[i] end

      local Qu = {}
      for i = 1, nu do Qu[i] = lu[i] + Bt_Vx[i] end

      -- Qxx = lxx + A^T * Vxx * A
      local At_Vxx = matrix.mul(At, Vxx)
      local At_Vxx_A = matrix.mul(At_Vxx, A)
      local Qxx = matrix.zeros(nx, nx)
      for i = 1, nx do
        for j = 1, nx do
          Qxx[i][j] = lxx[i][j] + At_Vxx_A[i][j]
        end
      end

      -- Quu = luu + B^T * Vxx * B + mu * I
      local Bt_Vxx = matrix.mul(Bt, Vxx)
      local Bt_Vxx_B = matrix.mul(Bt_Vxx, B)
      local Quu = matrix.zeros(nu, nu)
      for i = 1, nu do
        for j = 1, nu do
          Quu[i][j] = luu[i][j] + Bt_Vxx_B[i][j] + (i == j and mu or 0.0)
        end
      end

      -- Qux = lux + B^T * Vxx * A
      local Bt_Vxx_A = matrix.mul(Bt_Vxx, A)
      local Qux = matrix.zeros(nu, nx)
      for i = 1, nu do
        for j = 1, nx do
          Qux[i][j] = lux[i][j] + Bt_Vxx_A[i][j]
        end
      end

      -- Solve for feedforward k_ff = -Quu^-1 * Qu and feedback K = -Quu^-1 * Qux
      local neg_Qu = {}
      for i = 1, nu do neg_Qu[i] = -Qu[i] end

      local success_k, k_ff = pcall(function() return matrix.solve_lu(Quu, neg_Qu) end)
      if not success_k then
        backward_success = false
        break
      end

      local K_mat = matrix.zeros(nu, nx)
      for j = 1, nx do
        local col_Qux = {}
        for i = 1, nu do col_Qux[i] = -Qux[i][j] end
        local ok, col_K = pcall(function() return matrix.solve_lu(Quu, col_Qux) end)
        if not ok then
          backward_success = false
          break
        end
        for i = 1, nu do
          K_mat[i][j] = col_K[i]
        end
      end

      if not backward_success then
        break
      end

      k_traj[k + 1] = k_ff
      K_traj[k + 1] = K_mat

      -- Update Vx = Qx + K^T * Quu * k + K^T * Qu + Qux^T * k
      -- Update Vxx = Qxx + K^T * Quu * K + K^T * Qux + Qux^T * K
      local Kt = matrix.transpose(K_mat)
      local Quu_k = matrix.mul_vec(Quu, k_ff)
      local Kt_Qu = matrix.mul_vec(Kt, Qu)
      local Kt_Quu_k = matrix.mul_vec(Kt, Quu_k)
      local Quxt = matrix.transpose(Qux)
      local Quxt_k = matrix.mul_vec(Quxt, k_ff)

      local new_Vx = {}
      for i = 1, nx do
        new_Vx[i] = Qx[i] + Kt_Qu[i] + Kt_Quu_k[i] + Quxt_k[i]
      end
      Vx = new_Vx

      local Kt_Quu = matrix.mul(Kt, Quu)
      local Kt_Quu_K = matrix.mul(Kt_Quu, K_mat)
      local Kt_Qux = matrix.mul(Kt, Qux)
      local Quxt_K = matrix.mul(Quxt, K_mat)

      local new_Vxx = matrix.zeros(nx, nx)
      for i = 1, nx do
        for j = 1, nx do
          new_Vxx[i][j] = Qxx[i][j] + Kt_Quu_K[i][j] + Kt_Qux[i][j] + Quxt_K[i][j]
        end
      end
      Vxx = new_Vxx

      local k_Qu = 0.0
      for i = 1, nu do k_Qu = k_Qu + k_ff[i] * Qu[i] end
      local k_Quu_k = 0.0
      for i = 1, nu do k_Quu_k = k_Quu_k + k_ff[i] * Quu_k[i] end
      expected_cost_reduction = expected_cost_reduction - (k_Qu + 0.5 * k_Quu_k)
    end

    if not backward_success then
      mu = math.min(mu * delta_mu, mu_max)
    else
      -- Forward Pass with line search
      local step_accepted = false
      local best_x_traj = nil
      local best_u_traj = nil
      local best_cost = current_cost

      for _, alpha in ipairs(alphas) do
        local new_x_traj = { x0 }
        local new_u_traj = {}
        local line_cost = 0.0

        for k = 0, N - 1 do
          local x_curr = new_x_traj[k + 1]
          local x_nom = x_traj[k + 1]
          local u_nom = u_traj[k + 1]
          local k_ff = k_traj[k + 1]
          local K_mat = K_traj[k + 1]

          local dx = {}
          for i = 1, nx do dx[i] = x_curr[i] - x_nom[i] end
          local K_dx = matrix.mul_vec(K_mat, dx)

          local u_new = {}
          for i = 1, nu do
            u_new[i] = u_nom[i] + alpha * k_ff[i] + K_dx[i]
          end
          table.insert(new_u_traj, u_new)

          line_cost = line_cost + prob.cost(x_curr, u_new, k)
          local x_next = prob.dynamics(x_curr, u_new, k)
          table.insert(new_x_traj, x_next)
        end

        line_cost = line_cost + prob.terminal_cost(new_x_traj[N + 1])

        if line_cost < current_cost then
          step_accepted = true
          best_x_traj = new_x_traj
          best_u_traj = new_u_traj
          best_cost = line_cost
          break
        end
      end

      if step_accepted then
        local delta_J = current_cost - best_cost
        current_cost = best_cost
        x_traj = best_x_traj
        u_traj = best_u_traj
        table.insert(cost_history, current_cost)

        mu = math.max(mu / delta_mu, mu_min)

        if math.abs(delta_J) < tol then
          converged = true
          break
        end
      else
        mu = math.min(mu * delta_mu, mu_max)
        if mu >= mu_max then
          break
        end
      end
    end
  end

  return {
    x = x_traj,
    u = u_traj,
    K = K_traj,
    k = k_traj,
    cost = current_cost,
    cost_history = cost_history,
    converged = converged,
    iterations = iter_count
  }
end

return ilqr
