local aba = require("fulcrum.dynamics.aba")

local gl = {}

-- Gauss-Legendre 4th-order (2-stage implicit Runge-Kutta)
local sqrt3 = math.sqrt(3.0)
local c1 = 0.5 - sqrt3 / 6.0
local c2 = 0.5 + sqrt3 / 6.0

local a11 = 0.25
local a12 = 0.25 - sqrt3 / 6.0
local a21 = 0.25 + sqrt3 / 6.0
local a22 = 0.25

local b1 = 0.5
local b2 = 0.5

function gl.step4(m, q, q_dot, tau, dt, max_iter, tol, f_ext)
  -- 2-stage Gauss-Legendre implicit Runge-Kutta:
  -- Symplectic, A-stable, order 4
  max_iter = max_iter or 20
  tol = tol or 1e-11

  local dof = m.dof

  -- Initial stage guesses
  local a0 = aba.forward_dynamics(m, q, q_dot, tau, f_ext)
  local k1_q = {}
  local k1_v = {}
  local k2_q = {}
  local k2_v = {}

  for i = 1, dof do
    k1_q[i] = q_dot[i] or 0.0
    k1_v[i] = a0[i] or 0.0
    k2_q[i] = q_dot[i] or 0.0
    k2_v[i] = a0[i] or 0.0
  end

  -- Fixed-point iteration on stages
  for _ = 1, max_iter do
    local q1 = {}
    local v1 = {}
    local q2 = {}
    local v2 = {}

    for i = 1, dof do
      q1[i] = (q[i] or 0.0) + dt * (a11 * k1_q[i] + a12 * k2_q[i])
      v1[i] = (q_dot[i] or 0.0) + dt * (a11 * k1_v[i] + a12 * k2_v[i])

      q2[i] = (q[i] or 0.0) + dt * (a21 * k1_q[i] + a22 * k2_q[i])
      v2[i] = (q_dot[i] or 0.0) + dt * (a21 * k1_v[i] + a22 * k2_v[i])
    end

    local a1 = aba.forward_dynamics(m, q1, v1, tau, f_ext)
    local a2 = aba.forward_dynamics(m, q2, v2, tau, f_ext)

    local err = 0.0
    for i = 1, dof do
      local new_k1_q = v1[i]
      local new_k1_v = a1[i]
      local new_k2_q = v2[i]
      local new_k2_v = a2[i]

      err = err + math.abs(new_k1_q - k1_q[i]) + math.abs(new_k1_v - k1_v[i])
                + math.abs(new_k2_q - k2_q[i]) + math.abs(new_k2_v - k2_v[i])

      k1_q[i] = new_k1_q
      k1_v[i] = new_k1_v
      k2_q[i] = new_k2_q
      k2_v[i] = new_k2_v
    end

    if err < tol then break end
  end

  -- Final update: y_{k+1} = y_k + dt * (b1 * k1 + b2 * k2)
  local q_next = {}
  local qd_next = {}

  for i = 1, dof do
    q_next[i] = (q[i] or 0.0) + dt * (b1 * k1_q[i] + b2 * k2_q[i])
    qd_next[i] = (q_dot[i] or 0.0) + dt * (b1 * k1_v[i] + b2 * k2_v[i])
  end

  local a_next = aba.forward_dynamics(m, q_next, qd_next, tau, f_ext)
  return q_next, qd_next, a_next
end

return gl
