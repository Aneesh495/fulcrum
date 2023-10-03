-- examples/02_double_pendulum.lua
-- Demonstrates chaotic double pendulum dynamics integrated with the
-- 4th-order Forest-Ruth symplectic integrator and Featherstone ABA algorithm.
-- Logs high-frequency telemetry and exports trajectory to CSV.

local vec3 = require("fulcrum.math.vec3")
local matrix = require("fulcrum.math.matrix")
local pendulums = require("fulcrum.models.pendulums")
local crba = require("fulcrum.dynamics.crba")
local forest_ruth = require("fulcrum.integrators.forest_ruth")
local csv_logger = require("fulcrum.viz.csv_logger")
local terminal = require("fulcrum.viz.terminal")

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: CHAOTIC DOUBLE PENDULUM DYNAMICS ===" .. terminal.RESET)
print(terminal.rgb_fg(139, 148, 158) .. "4th-Order Symplectic Forest-Ruth Integration & Featherstone Articulated Dynamics\n" .. terminal.RESET)

-- Build 2-DOF double pendulum
local model, meta = pendulums.double_pendulum({
  m1 = 1.0,
  m2 = 1.0,
  l1 = 1.0,
  l2 = 1.0,
  gravity = vec3.new(0, 0, -9.81)
})

-- High-energy chaotic initial condition
local q = { 2.1, -1.6 }  -- rad
local qd = { 0.0, 0.0 }  -- rad/s
local tau = { 0.0, 0.0 } -- passive unforced

local kinematics = require("fulcrum.dynamics.kinematics")

local g = 9.81
local function compute_total_energy(q_curr, qd_curr)
  -- Kinetic energy: 0.5 * q_dot^T * M(q) * q_dot
  local M = crba.inertia_matrix(model, q_curr)
  local M_qd = matrix.mul_vec(M, qd_curr)
  local T = 0.5 * (qd_curr[1] * M_qd[1] + qd_curr[2] * M_qd[2])

  -- Potential energy: -mass * dot(gravity, com_world)
  local grav_3d = vec3.new(model.gravity[4], model.gravity[5], model.gravity[6])
  local com_w, total_mass = kinematics.center_of_mass(model, q_curr)
  local V = -total_mass * vec3.dot(grav_3d, com_w)

  return T + V
end

local E0 = compute_total_energy(q, qd)

local logger = csv_logger.new({
  "time", "theta1", "theta2", "omega1", "omega2", "energy", "energy_drift"
})

local dt = 0.002
local total_steps = 1000
local print_interval = 100

local max_drift = 0.0

for step = 1, total_steps do
  local t = step * dt
  q, qd = forest_ruth.step(model, q, qd, tau, dt)

  local E = compute_total_energy(q, qd)
  local drift = math.abs((E - E0) / E0)
  if drift > max_drift then max_drift = drift end

  logger:log({
    t, q[1], q[2], qd[1], qd[2], E, drift
  })

  if step % print_interval == 0 then
    local prog = step / total_steps
    local bar = terminal.progress_bar(prog, 20)
    io.write(string.format("\r[%s] Step %4d/%d (t=%.2fs) | q1=%+.3f q2=%+.3f | Max Drift: %.2e",
      bar, step, total_steps, t, q[1], q[2], max_drift))
    io.flush()
  end
end

print("\n")

local csv_path = "examples/double_pendulum_trajectory.csv"
logger:write_file(csv_path)

local summary = {
  string.format("Simulation Time:     %.3f seconds (%d steps at dt = %.4f s)", total_steps * dt, total_steps, dt),
  string.format("Initial Energy (E0): %+.6f J", E0),
  string.format("Final Energy:        %+.6f J", compute_total_energy(q, qd)),
  string.format("Max Relative Drift:  %.3e (Forest-Ruth 4th order symplectic)", max_drift),
  string.format("Exported Datafile:   %s", csv_path)
}
print(terminal.draw_box("SIMULATION SUMMARY", summary, 74))
print(terminal.rgb_fg(0, 230, 118) .. terminal.BOLD .. "Double pendulum simulation completed successfully." .. terminal.RESET .. "\n")
