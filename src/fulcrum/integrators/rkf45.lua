local aba = require("fulcrum.dynamics.aba")

local rkf45 = {}

function rkf45.step(m, q, q_dot, tau, dt, tol, f_ext)
  -- Embedded Runge-Kutta-Fehlberg 4(5) with adaptive step size
  tol = tol or 1e-6
  local dof = m.dof

  -- Stage evaluations for y = [q, v]
  -- dy/dt = [v, a(q, v)]
  local k1_q = {}
  local k1_v = aba.forward_dynamics(m, q, q_dot, tau, f_ext)
  for i = 1, dof do k1_q[i] = q_dot[i] or 0.0 end

  -- Stage 2: c2 = 1/4
  local q2 = {}
  local v2 = {}
  for i = 1, dof do
    q2[i] = (q[i] or 0.0) + dt * (0.25 * k1_q[i])
    v2[i] = (q_dot[i] or 0.0) + dt * (0.25 * k1_v[i])
  end
  local k2_v = aba.forward_dynamics(m, q2, v2, tau, f_ext)
  local k2_q = v2

  -- Stage 3: c3 = 3/8
  local q3 = {}
  local v3 = {}
  for i = 1, dof do
    q3[i] = (q[i] or 0.0) + dt * (0.09375 * k1_q[i] + 0.28125 * k2_q[i])
    v3[i] = (q_dot[i] or 0.0) + dt * (0.09375 * k1_v[i] + 0.28125 * k2_v[i])
  end
  local k3_v = aba.forward_dynamics(m, q3, v3, tau, f_ext)
  local k3_q = v3

  -- Stage 4: c4 = 12/13
  local q4 = {}
  local v4 = {}
  for i = 1, dof do
    q4[i] = (q[i] or 0.0) + dt * ((1932.0/2197.0) * k1_q[i] - (7200.0/2197.0) * k2_q[i] + (7296.0/2197.0) * k3_q[i])
    v4[i] = (q_dot[i] or 0.0) + dt * ((1932.0/2197.0) * k1_v[i] - (7200.0/2197.0) * k2_v[i] + (7296.0/2197.0) * k3_v[i])
  end
  local k4_v = aba.forward_dynamics(m, q4, v4, tau, f_ext)
  local k4_q = v4

  -- Stage 5: c5 = 1
  local q5 = {}
  local v5 = {}
  for i = 1, dof do
    q5[i] = (q[i] or 0.0) + dt * ((439.0/216.0) * k1_q[i] - 8.0 * k2_q[i] + (3680.0/513.0) * k3_q[i] - (845.0/4104.0) * k4_q[i])
    v5[i] = (q_dot[i] or 0.0) + dt * ((439.0/216.0) * k1_v[i] - 8.0 * k1_v[i] + (3680.0/513.0) * k3_v[i] - (845.0/4104.0) * k4_v[i])
  end
  local k5_v = aba.forward_dynamics(m, q5, v5, tau, f_ext)
  local k5_q = v5

  -- 4th order solution (y4) and 5th order solution (y5)
  local q_next = {}
  local v_next = {}
  local max_err = 0.0

  for i = 1, dof do
    -- 4th order update: b4 = [ 25/216, 0, 1408/2565, 2197/4104, -1/5 ]
    local q4_val = (q[i] or 0.0) + dt * ((25.0/216.0) * k1_q[i] + (1408.0/2565.0) * k3_q[i] + (2197.0/4104.0) * k4_q[i] - 0.2 * k5_q[i])
    local v4_val = (q_dot[i] or 0.0) + dt * ((25.0/216.0) * k1_v[i] + (1408.0/2565.0) * k3_v[i] + (2197.0/4104.0) * k4_v[i] - 0.2 * k5_v[i])

    -- 5th order update: b5 = [ 16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55 ]
    -- Error estimation: difference between 4th and 5th order
    local q5_val = q4_val -- approximation of difference
    local err_q = math.abs(dt * ((1.0/360.0) * k1_q[i] - (128.0/4275.0) * k3_q[i] - (2197.0/75240.0) * k4_q[i] + (1.0/50.0) * k5_q[i]))
    local err_v = math.abs(dt * ((1.0/360.0) * k1_v[i] - (128.0/4275.0) * k3_v[i] - (2197.0/75240.0) * k4_v[i] + (1.0/50.0) * k5_v[i]))

    max_err = math.max(max_err, err_q + err_v)
    q_next[i] = q4_val
    v_next[i] = v4_val
  end

  -- Optimal new step size
  local factor = 0.9 * ((tol / (max_err + 1e-15)) ^ 0.2)
  local dt_next = dt * math.max(0.2, math.min(2.0, factor))

  return q_next, v_next, dt_next, max_err
end

return rkf45
