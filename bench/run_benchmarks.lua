local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local builder = require("fulcrum.models.builder")
local models = require("fulcrum.models")
local aba = require("fulcrum.dynamics.aba")
local crba = require("fulcrum.dynamics.crba")
local rnea = require("fulcrum.dynamics.rnea")
local symplectic = require("fulcrum.integrators.symplectic")
local forest_ruth = require("fulcrum.integrators.forest_ruth")
local lgvi = require("fulcrum.integrators.lie_group_integrator")
local lemke = require("fulcrum.contact.lemke")
local pgs = require("fulcrum.contact.pgs")
local terminal = require("fulcrum.viz.terminal")

local function print_banner(title)
  local width = 76
  print("\n" .. terminal.rgb_fg(0, 229, 255) .. terminal.BOLD)
  print(string.rep("═", width))
  print("  " .. title)
  print(string.rep("═", width) .. terminal.RESET)
end

local function print_header(cols)
  local parts = {}
  for _, col in ipairs(cols) do
    table.insert(parts, string.format("%-" .. col.width .. "s", col.name))
  end
  print(terminal.rgb_fg(139, 148, 158) .. table.concat(parts, " | ") .. terminal.RESET)
  print(string.rep("─", 76))
end

local function print_row(vals, widths)
  local parts = {}
  for i, val in ipairs(vals) do
    table.insert(parts, string.format("%-" .. widths[i] .. "s", tostring(val)))
  end
  print(table.concat(parts, " | "))
end

-- Benchmark 1: Forward Dynamics Scaling (ABA vs CRBA + LU Solve)
local function bench_forward_dynamics()
  print_banner("Benchmark 1: Featherstone O(N) ABA vs O(N^3) CRBA + LU Inversion")

  local dof_list = { 2, 4, 7, 12, 18, 25 }
  local widths = { 8, 14, 16, 14, 12 }

  print_header({
    { name = "DOF (N)", width = widths[1] },
    { name = "ABA (us/step)", width = widths[2] },
    { name = "CRBA+LU (us/step)", width = widths[3] },
    { name = "Speedup", width = widths[4] },
    { name = "Agreement", width = widths[5] }
  })

  for _, n_dof in ipairs(dof_list) do
    -- Create serial chain of n_dof revolute joints
    local b = builder.new("chain_" .. n_dof, vec3.new(0, 0, -9.81))
    b:add_base("world", { mass = 0.0 })
    for i = 1, n_dof do
      local p_name = (i == 1) and "world" or ("link_" .. (i - 1))
      local axis = (i % 2 == 1) and vec3.unit_z() or vec3.unit_y()
      b:add_link({
        name = "link_" .. i,
        parent = p_name,
        joint_type = "revolute",
        axis = axis,
        origin = { xyz = vec3.new(0, 0, 0.2), rpy = vec3.zero() },
        mass = 1.0,
        com = vec3.new(0, 0, 0.1),
        inertia = { ixx = 0.02, iyy = 0.02, izz = 0.005 }
      })
    end
    local model = b:build()

    local q = {}
    local qd = {}
    local tau = {}
    for i = 1, n_dof do
      q[i] = 0.2 * math.sin(i * 0.5)
      qd[i] = 0.1 * math.cos(i * 0.5)
      tau[i] = 1.0 * math.sin(i * 0.3)
    end

    -- Warmup
    local qdd_aba = aba.forward_dynamics(model, q, qd, tau)

    -- Benchmark ABA
    local iters = 100
    local t0 = os.clock()
    for _ = 1, iters do
      aba.forward_dynamics(model, q, qd, tau)
    end
    local time_aba_us = ((os.clock() - t0) / iters) * 1e6

    -- Benchmark CRBA + RNEA bias + LU Solve
    local t1 = os.clock()
    local qdd_crba
    for _ = 1, iters do
      local M = crba.inertia_matrix(model, q)
      local qdd_zero = {}
      for i = 1, n_dof do qdd_zero[i] = 0.0 end
      local b_bias = rnea.inverse_dynamics(model, q, qd, qdd_zero)
      local rhs = {}
      for i = 1, n_dof do rhs[i] = tau[i] - b_bias[i] end
      qdd_crba = matrix.solve_lu(M, rhs)
    end
    local time_crba_us = ((os.clock() - t1) / iters) * 1e6

    local max_err = 0.0
    for i = 1, n_dof do
      local err = math.abs(qdd_aba[i] - qdd_crba[i])
      if err > max_err then max_err = err end
    end

    local speedup = string.format("%.2fx", time_crba_us / time_aba_us)
    local agree = (max_err < 1e-5) and "Exact (1e-6)" or "Fail"

    print_row({
      tostring(n_dof),
      string.format("%.1f", time_aba_us),
      string.format("%.1f", time_crba_us),
      speedup,
      agree
    }, widths)
  end
end

-- Benchmark 2: Symplectic Energy Conservation over 10,000 steps
local function bench_symplectic_conservation()
  print_banner("Benchmark 2: Symplectic vs Classical Integrator Energy Drift (10,000 Steps)")

  local dt = 0.002
  local steps = 10000

  local widths = { 26, 12, 16, 16 }
  print_header({
    { name = "Integrator Method", width = widths[1] },
    { name = "Order", width = widths[2] },
    { name = "Energy Drift", width = widths[3] },
    { name = "Runtime (ms)", width = widths[4] }
  })

  -- Hamiltonian system: nonlinear pendulum dth/dt = w, dw/dt = -sin(th)
  local function pendulum_accel(th) return -math.sin(th) end
  local function pendulum_energy(th, w) return 0.5 * w * w - math.cos(th) end

  local th0 = 1.8 -- high-amplitude nonlinear swing
  local w0 = 0.0
  local E0 = pendulum_energy(th0, w0)

  -- 1. Standard Explicit Euler
  local th, w = th0, w0
  local t0 = os.clock()
  local max_drift_euler = 0.0
  for _ = 1, steps do
    local a = pendulum_accel(th)
    th = th + w * dt
    w = w + a * dt
    local d = math.abs((pendulum_energy(th, w) - E0) / E0)
    if d > max_drift_euler then max_drift_euler = d end
  end
  local time_euler_ms = (os.clock() - t0) * 1e3
  print_row({ "Explicit Euler (Forward)", "1st", string.format("%.2e", max_drift_euler), string.format("%.1f", time_euler_ms) }, widths)

  -- 2. Classical Runge-Kutta 4th-Order (RK4)
  th, w = th0, w0
  t0 = os.clock()
  local max_drift_rk4 = 0.0
  for _ = 1, steps do
    local k1_th = w
    local k1_w = pendulum_accel(th)
    local k2_th = w + 0.5 * dt * k1_w
    local k2_w = pendulum_accel(th + 0.5 * dt * k1_th)
    local k3_th = w + 0.5 * dt * k2_w
    local k3_w = pendulum_accel(th + 0.5 * dt * k2_th)
    local k4_th = w + dt * k3_w
    local k4_w = pendulum_accel(th + dt * k3_th)

    th = th + (dt / 6.0) * (k1_th + 2 * k2_th + 2 * k3_th + k4_th)
    w = w + (dt / 6.0) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)

    local d = math.abs((pendulum_energy(th, w) - E0) / E0)
    if d > max_drift_rk4 then max_drift_rk4 = d end
  end
  local time_rk4_ms = (os.clock() - t0) * 1e3
  print_row({ "Classical Runge-Kutta 4", "4th", string.format("%.2e", max_drift_rk4), string.format("%.1f", time_rk4_ms) }, widths)

  -- 3. Symplectic Velocity Verlet
  th, w = th0, w0
  t0 = os.clock()
  local max_drift_verlet = 0.0
  for _ = 1, steps do
    local q_next, qd_next = symplectic.velocity_verlet_step({ th }, { w }, function(q_in) return { pendulum_accel(q_in[1]) } end, dt)
    th, w = q_next[1], qd_next[1]
    local d = math.abs((pendulum_energy(th, w) - E0) / E0)
    if d > max_drift_verlet then max_drift_verlet = d end
  end
  local time_verlet_ms = (os.clock() - t0) * 1e3
  print_row({ "Symplectic Velocity Verlet", "2nd", string.format("%.2e", max_drift_verlet), string.format("%.1f", time_verlet_ms) }, widths)

  -- 4. Forest-Ruth 4th-Order Symplectic
  th, w = th0, w0
  t0 = os.clock()
  local max_drift_fr = 0.0
  for _ = 1, steps do
    local q_next, qd_next = forest_ruth.step_hamiltonian({ th }, { w }, function(q_in) return { pendulum_accel(q_in[1]) } end, dt)
    th, w = q_next[1], qd_next[1]
    local d = math.abs((pendulum_energy(th, w) - E0) / E0)
    if d > max_drift_fr then max_drift_fr = d end
  end
  local time_fr_ms = (os.clock() - t0) * 1e3
  print_row({ "Forest-Ruth Symplectic", "4th", string.format("%.2e", max_drift_fr), string.format("%.1f", time_fr_ms) }, widths)

  -- 5. Lie Group Variational Integrator on SO(3) Free Flight
  local J_rigid = vec3.new(1.0, 2.0, 3.0)
  local w_rigid = vec3.new(0.01, 3.0, 0.0)
  local R_rigid = mat3.identity()
  local E0_rigid = 0.5 * (1.0 * w_rigid.x^2 + 2.0 * w_rigid.y^2 + 3.0 * w_rigid.z^2)

  t0 = os.clock()
  local max_drift_lgvi = 0.0
  for _ = 1, steps do
    R_rigid, w_rigid = lgvi.step_so3(R_rigid, w_rigid, J_rigid, nil, dt)
    local E_cur = 0.5 * (1.0 * w_rigid.x^2 + 2.0 * w_rigid.y^2 + 3.0 * w_rigid.z^2)
    local d = math.abs((E_cur - E0_rigid) / E0_rigid)
    if d > max_drift_lgvi then max_drift_lgvi = d end
  end
  local time_lgvi_ms = (os.clock() - t0) * 1e3
  print_row({ "Lie Group Variational (SO3)", "Symplectic", string.format("%.2e", max_drift_lgvi), string.format("%.1f", time_lgvi_ms) }, widths)
end

-- Benchmark 3: Contact LCP Solvers (Lemke vs Projected Gauss-Seidel)
local function bench_contact_solvers()
  print_banner("Benchmark 3: Contact Complementarity Solvers (Lemke vs Projected Gauss-Seidel)")

  local contact_counts = { 2, 4, 8, 16 }
  local widths = { 12, 14, 16, 14, 16 }

  print_header({
    { name = "Contacts (k)", width = widths[1] },
    { name = "Lemke (us)", width = widths[2] },
    { name = "Lemke Residual", width = widths[3] },
    { name = "PGS (us)", width = widths[4] },
    { name = "PGS Residual", width = widths[5] }
  })

  for _, nc in ipairs(contact_counts) do
    -- Construct positive-definite Delassus contact matrix
    local M = matrix.zeros(nc, nc)
    local q_vec = {}
    for i = 1, nc do
      q_vec[i] = -1.5 - 0.5 * math.sin(i)
      for j = 1, nc do
        if i == j then
          M[i][j] = 2.0 + (i * 0.1)
        else
          M[i][j] = 0.3 / (1.0 + math.abs(i - j))
        end
      end
    end

    local iters = 200

    -- Benchmark Lemke
    local t0 = os.clock()
    local z_lemke
    for _ = 1, iters do
      local _, sol = lemke.solve(M, q_vec)
      z_lemke = sol
    end
    local time_lemke_us = ((os.clock() - t0) / iters) * 1e6

    -- Compute Lemke residual: w = M*z + q, check w^T * z
    local comp_lemke = 0.0
    if z_lemke then
      local M_z = matrix.mul_vec(M, z_lemke)
      for i = 1, nc do
        local wi = M_z[i] + q_vec[i]
        comp_lemke = comp_lemke + math.abs(wi * z_lemke[i])
      end
    end

    -- Benchmark PGS
    local t1 = os.clock()
    local z_pgs
    for _ = 1, iters do
      z_pgs = pgs.solve_unilateral(M, q_vec, 100, 1e-5)
    end
    local time_pgs_us = ((os.clock() - t1) / iters) * 1e6

    local comp_pgs = 0.0
    if z_pgs then
      local M_z = matrix.mul_vec(M, z_pgs)
      for i = 1, nc do
        local wi = M_z[i] + q_vec[i]
        comp_pgs = comp_pgs + math.abs(wi * z_pgs[i])
      end
    end

    print_row({
      tostring(nc),
      string.format("%.1f", time_lemke_us),
      string.format("%.2e", comp_lemke),
      string.format("%.1f", time_pgs_us),
      string.format("%.2e", comp_pgs)
    }, widths)
  end
end

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM PERFORMANCE & PHYSICS BENCHMARK SUITE ===" .. terminal.RESET)
bench_forward_dynamics()
bench_symplectic_conservation()
bench_contact_solvers()
print("\n" .. terminal.rgb_fg(0, 230, 118) .. "All benchmarks executed successfully." .. terminal.RESET .. "\n")
