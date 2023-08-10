local friction = require("fulcrum.contact.friction")

local pgs = {}

function pgs.solve_box(M, q, lb, ub, max_iter, tol, omega, z_init)
  -- Projected Gauss-Seidel for box-constrained LCP: lb <= z <= ub
  max_iter = max_iter or 100
  tol = tol or 1e-6
  omega = omega or 1.2 -- SOR over-relaxation parameter
  local n = #q

  local z = {}
  for i = 1, n do
    z[i] = z_init and z_init[i] or 0.0
  end

  local iter = 0
  for step = 1, max_iter do
    iter = step
    local max_change = 0.0

    for i = 1, n do
      local M_ii = M[i][i]
      if math.abs(M_ii) > 1e-15 then
        local sum = q[i]
        for j = 1, n do
          if j ~= i then
            sum = sum + M[i][j] * z[j]
          end
        end

        local z_unconstrained = -sum / M_ii
        local z_sor = (1.0 - omega) * z[i] + omega * z_unconstrained

        -- Box projection
        local low = lb and lb[i] or 0.0
        local high = ub and ub[i] or 1e20
        local z_proj = math.max(low, math.min(high, z_sor))

        local diff = math.abs(z_proj - z[i])
        max_change = math.max(max_change, diff)
        z[i] = z_proj
      end
    end

    if max_change < tol then break end
  end

  return z, iter
end

function pgs.solve_contact_friction(W, b, num_contacts, mu, max_iter, tol, omega)
  -- Solves coupled 3D normal + tangential friction contact problem:
  -- For contact k, z has [fn, ft1, ft2].
  -- W * z + b is relative contact velocity.
  max_iter = max_iter or 100
  tol = tol or 1e-6
  omega = omega or 1.0
  mu = mu or 0.5

  local n = num_contacts * 3
  local z = {}
  for i = 1, n do z[i] = 0.0 end

  for step = 1, max_iter do
    local max_change = 0.0

    for c = 1, num_contacts do
      local idx_n  = (c - 1) * 3 + 1
      local idx_t1 = (c - 1) * 3 + 2
      local idx_t2 = (c - 1) * 3 + 3

      -- Normal update
      local M_nn = W[idx_n][idx_n]
      if math.abs(M_nn) > 1e-15 then
        local sum_n = b[idx_n]
        for j = 1, n do
          if j ~= idx_n then sum_n = sum_n + W[idx_n][j] * z[j] end
        end
        local zn_un = -sum_n / M_nn
        local zn_new = math.max(0.0, (1.0 - omega) * z[idx_n] + omega * zn_un)
        max_change = math.max(max_change, math.abs(zn_new - z[idx_n]))
        z[idx_n] = zn_new
      end

      -- Tangential updates
      local M_t1 = W[idx_t1][idx_t1]
      local zt1_new = z[idx_t1]
      if math.abs(M_t1) > 1e-15 then
        local sum_t1 = b[idx_t1]
        for j = 1, n do
          if j ~= idx_t1 then sum_t1 = sum_t1 + W[idx_t1][j] * z[j] end
        end
        local zt1_un = -sum_t1 / M_t1
        zt1_new = (1.0 - omega) * z[idx_t1] + omega * zt1_un
      end

      local M_t2 = W[idx_t2][idx_t2]
      local zt2_new = z[idx_t2]
      if math.abs(M_t2) > 1e-15 then
        local sum_t2 = b[idx_t2]
        for j = 1, n do
          if j ~= idx_t2 then sum_t2 = sum_t2 + W[idx_t2][j] * z[j] end
        end
        local zt2_un = -sum_t2 / M_t2
        zt2_new = (1.0 - omega) * z[idx_t2] + omega * zt2_un
      end

      -- Project onto circular Coulomb cone
      local fn, ft1, ft2 = friction.project_coulomb_cone(z[idx_n], zt1_new, zt2_new, mu)
      max_change = math.max(max_change, math.abs(fn - z[idx_n]))
      max_change = math.max(max_change, math.abs(ft1 - z[idx_t1]))
      max_change = math.max(max_change, math.abs(ft2 - z[idx_t2]))

      z[idx_n] = fn
      z[idx_t1] = ft1
      z[idx_t2] = ft2
    end

    if max_change < tol then break end
  end

  return z
end

return pgs
