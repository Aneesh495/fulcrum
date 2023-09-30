local aba = require("fulcrum.dynamics.aba")

local symplectic = {}

function symplectic.euler(m, q, q_dot, tau, dt, f_ext)
  -- Symplectic Euler (semi-implicit):
  -- v_{k+1} = v_k + dt * a(q_k, v_k)
  -- q_{k+1} = q_k + dt * v_{k+1}
  local qdd = aba.forward_dynamics(m, q, q_dot, tau, f_ext)

  local qd_next = {}
  local q_next = {}

  for i = 1, m.dof do
    local v_new = (q_dot[i] or 0.0) + (qdd[i] or 0.0) * dt
    qd_next[i] = v_new
    q_next[i] = (q[i] or 0.0) + v_new * dt
  end

  return q_next, qd_next, qdd
end

function symplectic.verlet(m, q, q_dot, tau, dt, f_ext)
  -- Velocity Verlet (Störmer-Verlet):
  -- 1. a_k = ABA(q_k, v_k)
  -- 2. v_{k+1/2} = v_k + 0.5 * dt * a_k
  -- 3. q_{k+1} = q_k + dt * v_{k+1/2}
  -- 4. a_{k+1} = ABA(q_{k+1}, v_{k+1/2})
  -- 5. v_{k+1} = v_{k+1/2} + 0.5 * dt * a_{k+1}
  local qdd_k = aba.forward_dynamics(m, q, q_dot, tau, f_ext)

  local qd_half = {}
  local q_next = {}

  for i = 1, m.dof do
    local v_half = (q_dot[i] or 0.0) + 0.5 * dt * qdd_k[i]
    qd_half[i] = v_half
    q_next[i] = (q[i] or 0.0) + dt * v_half
  end

  local qdd_next = aba.forward_dynamics(m, q_next, qd_half, tau, f_ext)

  local qd_next = {}
  for i = 1, m.dof do
    qd_next[i] = qd_half[i] + 0.5 * dt * qdd_next[i]
  end

  return q_next, qd_next, qdd_next
end

function symplectic.implicit_midpoint(m, q, q_dot, tau, dt, max_iter, tol, f_ext)
  -- Implicit Midpoint Rule (2nd order, symplectic, preserves all quadratic invariants)
  -- Solves:
  -- q_{k+1} = q_k + dt * 0.5 * (v_k + v_{k+1})
  -- v_{k+1} = v_k + dt * a(0.5*(q_k + q_{k+1}), 0.5*(v_k + v_{k+1}))
  max_iter = max_iter or 15
  tol = tol or 1e-10

  -- Initial guess from explicit Euler
  local qdd0 = aba.forward_dynamics(m, q, q_dot, tau, f_ext)
  local q1 = {}
  local v1 = {}
  for i = 1, m.dof do
    v1[i] = (q_dot[i] or 0.0) + dt * qdd0[i]
    q1[i] = (q[i] or 0.0) + dt * v1[i]
  end

  -- Fixed-point / Newton iterations
  for _ = 1, max_iter do
    local q_mid = {}
    local v_mid = {}
    for i = 1, m.dof do
      q_mid[i] = 0.5 * ((q[i] or 0.0) + q1[i])
      v_mid[i] = 0.5 * ((q_dot[i] or 0.0) + v1[i])
    end

    local a_mid = aba.forward_dynamics(m, q_mid, v_mid, tau, f_ext)

    local err = 0.0
    for i = 1, m.dof do
      local v1_new = (q_dot[i] or 0.0) + dt * a_mid[i]
      local q1_new = (q[i] or 0.0) + dt * 0.5 * ((q_dot[i] or 0.0) + v1_new)
      err = err + math.abs(v1_new - v1[i]) + math.abs(q1_new - q1[i])
      v1[i] = v1_new
      q1[i] = q1_new
    end

    if err < tol then break end
  end

  local a_final = aba.forward_dynamics(m, q1, v1, tau, f_ext)
  return q1, v1, a_final
end

function symplectic.verlet_step(q, qd, accel_fn, dt)
  local n = #q
  local a_k = accel_fn(q)
  local qd_half = {}
  local q_next = {}
  for i = 1, n do
    local v_half = (qd[i] or 0.0) + 0.5 * dt * (a_k[i] or 0.0)
    qd_half[i] = v_half
    q_next[i] = (q[i] or 0.0) + dt * v_half
  end

  local a_next = accel_fn(q_next)
  local qd_next = {}
  for i = 1, n do
    qd_next[i] = qd_half[i] + 0.5 * dt * (a_next[i] or 0.0)
  end

  return q_next, qd_next, a_next
end

symplectic.velocity_verlet_step = symplectic.verlet_step

return symplectic
