local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")
local rnea = require("fulcrum.dynamics.rnea")
local crba = require("fulcrum.dynamics.crba")

local bilateral = {}

function bilateral.solve_constrained_acceleration(m, q, q_dot, tau, G, g_val, g_dot_val, alpha, beta, f_ext)
  -- Solves KKT system with Baumgarte stabilization:
  -- [ M,   -G^T ] [ q_ddot ] = [ tau - C*qd - g(q) ]
  -- [ G,    0   ] [ lambda ]   [ -2*alpha*g_dot - beta^2*g ]
  alpha = alpha or 10.0 -- damping frequency
  beta = beta or 10.0   -- natural frequency

  local dof = m.dof
  local k = G.rows

  local M = crba.inertia_matrix(m, q)
  local c = rnea.coriolis_centrifugal(m, q, q_dot)
  local grav = rnea.gravity_terms(m, q)

  local tau_eff = {}
  for i = 1, dof do
    local applied = tau and tau[i] or 0.0
    tau_eff[i] = applied - c[i] - grav[i]
  end

  local rhs_constraint = {}
  for i = 1, k do
    local g_i = g_val and g_val[i] or 0.0
    local gd_i = g_dot_val and g_dot_val[i] or 0.0
    rhs_constraint[i] = -2.0 * alpha * gd_i - (beta * beta) * g_i
  end

  -- Assemble KKT matrix of size (dof + k) x (dof + k)
  local KKT = matrix.zeros(dof + k, dof + k)
  local rhs = {}

  for i = 1, dof do
    for j = 1, dof do
      KKT[i][j] = M[i][j]
    end
    for j = 1, k do
      KKT[i][dof + j] = G[j][i]
      KKT[dof + j][i] = G[j][i]
    end
    rhs[i] = tau_eff[i]
  end

  for j = 1, k do
    rhs[dof + j] = rhs_constraint[j]
  end

  -- Regularize lower-right constraint block slightly to ensure solvability
  for j = 1, k do
    KKT[dof + j][dof + j] = -1e-8
  end

  local sol = matrix.solve_lu(KKT, rhs)

  local q_ddot = {}
  local lambda = {}
  for i = 1, dof do
    q_ddot[i] = sol[i]
  end
  for j = 1, k do
    lambda[j] = sol[dof + j]
  end

  return q_ddot, lambda
end

return bilateral
