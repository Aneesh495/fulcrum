-- examples/07_trajectory_optimization.lua
-- Demonstrates non-linear optimal control and trajectory optimization using:
-- 1. Iterative Linear Quadratic Regulator (iLQR / DDP) with Riccati recursion
-- 2. Hermite-Simpson Direct Collocation transcribed into a Sequential Quadratic Program (SQP)

local ilqr = require("fulcrum.opt.ilqr")
local collocation = require("fulcrum.opt.collocation")
local terminal = require("fulcrum.viz.terminal")

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: OPTIMAL TRAJECTORY GENERATION ===" .. terminal.RESET)
print(terminal.rgb_fg(139, 148, 158) .. "Differential Dynamic Programming (iLQR) & Hermite-Simpson Direct Collocation\n" .. terminal.RESET)

-- =========================================================================
-- Part 1: iLQR Non-Linear Pendulum Stabilization
-- =========================================================================
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "--- Part 1: iLQR / DDP Dynamic Programming ---" .. terminal.RESET)

local dt_ilqr = 0.05
local N_ilqr = 30
local g = 9.81
local l = 1.0

local ilqr_problem = {
  n_x = 2,
  n_u = 1,
  N = N_ilqr,
  dt = dt_ilqr,
  dynamics = function(x, u, _)
    local th, w = x[1], x[2]
    local tau = u[1]
    local w_dot = -(g / l) * math.sin(th) + tau
    return {
      th + w * dt_ilqr,
      w + w_dot * dt_ilqr
    }
  end,
  cost = function(x, u, _)
    local th, w = x[1], x[2]
    local tau = u[1]
    return 0.5 * (25.0 * th * th + 2.5 * w * w + 0.1 * tau * tau)
  end,
  terminal_cost = function(x)
    local th, w = x[1], x[2]
    return 0.5 * (150.0 * th * th + 30.0 * w * w)
  end
}

local x0_ilqr = { 0.85, 0.0 } -- Displaced angle 0.85 rad (~49 deg)
local res_ilqr = ilqr.solve(ilqr_problem, x0_ilqr, nil, {
  max_iter = 35,
  tol = 1e-4
})

local x_final_ilqr = res_ilqr.x[N_ilqr + 1]
local cost_init = res_ilqr.cost_history[1]
local cost_final = res_ilqr.cost

print(string.format("  Initial State:       theta = %+.3f rad,  omega = %+.3f rad/s", x0_ilqr[1], x0_ilqr[2]))
print(string.format("  Final State:         theta = %+.4f rad, omega = %+.4f rad/s", x_final_ilqr[1], x_final_ilqr[2]))
print(string.format("  Initial Cost (J0):   %.4f", cost_init))
print(string.format("  Optimized Cost (J*): %.4f (Reduced by %.1f%%)", cost_final, (1.0 - cost_final / cost_init) * 100.0))
print(string.format("  Iterations:          %d Riccati backward/forward passes\n", #res_ilqr.cost_history - 1))

-- =========================================================================
-- Part 2: Hermite-Simpson Direct Collocation with SQP
-- =========================================================================
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "--- Part 2: Hermite-Simpson Direct Collocation ---" .. terminal.RESET)

local segments = 8
local duration = 1.0

local colloc_problem = {
  n_x = 2,
  n_u = 1,
  segments = segments,
  duration = duration,
  x_start = { 0.0, 0.0 },
  x_goal = { 1.0, 0.0 },
  dynamics = function(x, u)
    local v = x[2]
    local a = u[1]
    return { v, a }
  end,
  cost = function(_, u)
    local a = u[1]
    return 0.5 * a * a
  end
}

local res_colloc = collocation.solve(colloc_problem, {
  max_iter = 12,
  tol_defect = 1e-5
})

local x_final_colloc = res_colloc.x[#res_colloc.x]
local max_defect = 0.0
for k = 1, segments do
  local def = collocation.hermite_simpson_defect(
    res_colloc.x[k], res_colloc.u[k], res_colloc.x[k + 1], res_colloc.u[k + 1],
    colloc_problem.dynamics, res_colloc.dt
  )
  for _, d in ipairs(def) do
    if math.abs(d) > max_defect then max_defect = math.abs(d) end
  end
end

print(string.format("  Boundary Conditions: x(0) = [0.0, 0.0] -> x(1) = [%.4f, %.4f]", x_final_colloc[1], x_final_colloc[2]))
print(string.format("  Collocation Segments:%d (Hermite-Simpson 4th order transcription)", segments))
print(string.format("  Max Kinematic Defect:%.2e (Equating mid-point Simpson polynomials)", max_defect))
print(string.format("  Total Effort Cost:   %.4f (Analytical minimum: 6.0000)\n", res_colloc.cost))

local summary = {
  string.format("iLQR Final Angle Error:     %.2e rad (Target: 0.0)", math.abs(x_final_ilqr[1])),
  string.format("iLQR Final Velocity Error:  %.2e rad/s (Target: 0.0)", math.abs(x_final_ilqr[2])),
  string.format("Collocation Final Pos Error:%.2e m (Target: 1.0)", math.abs(x_final_colloc[1] - 1.0)),
  string.format("Collocation Max Defect:     %.2e (KKT constraint violation)", max_defect)
}
print(terminal.draw_box("OPTIMAL CONTROL SUMMARY", summary, 76))
print(terminal.rgb_fg(0, 230, 118) .. terminal.BOLD .. "Optimal trajectory generation demonstrations completed successfully." .. terminal.RESET .. "\n")
