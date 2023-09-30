local aba = require("fulcrum.dynamics.aba")

local forest_ruth = {}

-- Forest-Ruth 4th-order symplectic coefficients
local cbrt2 = 2.0 ^ (1.0 / 3.0)
local theta = 1.0 / (2.0 - cbrt2)

local c = {
  0.5 * theta,
  0.5 * (1.0 - theta),
  0.5 * (1.0 - theta),
  0.5 * theta
}

local d = {
  theta,
  1.0 - 2.0 * theta,
  theta,
  0.0
}

function forest_ruth.step(m, q, q_dot, tau, dt, f_ext)
  -- 4th-order symplectic integration step
  local q_curr = {}
  local qd_curr = {}

  for i = 1, m.dof do
    q_curr[i] = q[i] or 0.0
    qd_curr[i] = q_dot[i] or 0.0
  end

  local qdd_last = nil

  for s = 1, 4 do
    -- Drift step: q = q + c_s * dt * qd
    if c[s] ~= 0.0 then
      for i = 1, m.dof do
        q_curr[i] = q_curr[i] + c[s] * dt * qd_curr[i]
      end
    end

    -- Kick step: qd = qd + d_s * dt * a
    if d[s] ~= 0.0 then
      qdd_last = aba.forward_dynamics(m, q_curr, qd_curr, tau, f_ext)
      for i = 1, m.dof do
        qd_curr[i] = qd_curr[i] + d[s] * dt * qdd_last[i]
      end
    end
  end

  qdd_last = aba.forward_dynamics(m, q_curr, qd_curr, tau, f_ext)
  return q_curr, qd_curr, qdd_last
end

-- Yoshida 6th-order symplectic integrator composition
local w1 = 1.0 / (2.0 - 2.0 ^ (1.0 / 5.0))
local w0 = 1.0 - 2.0 * w1

function forest_ruth.yoshida6(m, q, q_dot, tau, dt, f_ext)
  local q1, v1 = forest_ruth.step(m, q, q_dot, tau, w1 * dt, f_ext)
  local q2, v2 = forest_ruth.step(m, q1, v1, tau, w0 * dt, f_ext)
  local q3, v3, a3 = forest_ruth.step(m, q2, v2, tau, w1 * dt, f_ext)
  return q3, v3, a3
end

function forest_ruth.step_hamiltonian(q, q_dot, accel_fn, dt)
  local n = #q
  local q_curr = {}
  local qd_curr = {}
  for i = 1, n do
    q_curr[i] = q[i] or 0.0
    qd_curr[i] = q_dot[i] or 0.0
  end

  for s = 1, 4 do
    if c[s] ~= 0.0 then
      for i = 1, n do
        q_curr[i] = q_curr[i] + c[s] * dt * qd_curr[i]
      end
    end
    if d[s] ~= 0.0 then
      local acc = accel_fn(q_curr)
      for i = 1, n do
        qd_curr[i] = qd_curr[i] + d[s] * dt * (acc[i] or 0.0)
      end
    end
  end

  return q_curr, qd_curr
end

return forest_ruth
