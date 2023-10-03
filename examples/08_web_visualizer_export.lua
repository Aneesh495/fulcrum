-- examples/08_web_visualizer_export.lua
-- Generates a multi-link articulated simulation trajectory with forward kinematics
-- and exports structured JSON telemetry to web/trajectory.json for browser playback.

local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local pendulums = require("fulcrum.models.pendulums")
local crba = require("fulcrum.dynamics.crba")
local kinematics = require("fulcrum.dynamics.kinematics")
local forest_ruth = require("fulcrum.integrators.forest_ruth")
local json_exporter = require("fulcrum.viz.json_exporter")
local terminal = require("fulcrum.viz.terminal")

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: WEB VISUALIZER TRAJECTORY EXPORT ===" .. terminal.RESET)
print(terminal.rgb_fg(139, 148, 158) .. "Recording Multi-Link Articulated Simulation to web/trajectory.json\n" .. terminal.RESET)

-- Double pendulum physical system
local model, meta = pendulums.double_pendulum({
  m1 = 1.0,
  m2 = 0.8,
  l1 = 0.8,
  l2 = 0.6,
  gravity = vec3.new(0.0, 0.0, -9.81)
})

local q = { 2.2, -1.8 }
local qd = { 0.0, 0.0 }
local tau = { 0.0, 0.0 }

local dt = 0.005
local total_steps = 400
local record_stride = 2 -- Record every 2nd step (100 Hz output)

local m1, m2 = meta.m1, meta.m2
local l1, l2 = meta.l1, meta.l2
local c1, c2 = 0.5 * l1, 0.5 * l2
local g = 9.81

local function compute_energy(q_curr, qd_curr)
  local M = crba.inertia_matrix(model, q_curr)
  local M_qd = matrix.mul_vec(M, qd_curr)
  local T = 0.5 * (qd_curr[1] * M_qd[1] + qd_curr[2] * M_qd[2])
  local V = m1 * g * c1 * math.cos(q_curr[1]) + m2 * g * (l1 * math.cos(q_curr[1]) + c2 * math.cos(q_curr[1] + q_curr[2]))
  return T + V, T, V
end

local E0, T0, V0 = compute_energy(q, qd)

local frames = {}

for step = 1, total_steps do
  local t = step * dt
  q, qd = forest_ruth.step(model, q, qd, tau, dt)

  if step % record_stride == 0 then
    local E, T, V = compute_energy(q, qd)
    local drift = math.abs((E - E0) / (math.abs(E0) > 1e-6 and math.abs(E0) or 1.0))

    -- Forward kinematics for link transforms
    local _, X_world = kinematics.forward_kinematics(model, q)
    local link_data = {}

    for b_id = 1, model.num_bodies do
      local b = model.bodies[b_id]
      local X_inv = matrix.inv(X_world[b_id])
      local R = mat3.new(
        X_inv[1][1], X_inv[1][2], X_inv[1][3],
        X_inv[2][1], X_inv[2][2], X_inv[2][3],
        X_inv[3][1], X_inv[3][2], X_inv[3][3]
      )
      table.insert(link_data, {
        name = b.name,
        pos = { -X_inv[5][3], X_inv[4][3], -X_inv[4][2] },
        rot = {
          { R[1][1], R[1][2], R[1][3] },
          { R[2][1], R[2][2], R[2][3] },
          { R[3][1], R[3][2], R[3][3] }
        }
      })
    end

    table.insert(frames, {
      time = t,
      q = { q[1], q[2] },
      qd = { qd[1], qd[2] },
      energy = E,
      kinetic = T,
      potential = V,
      drift = drift,
      links = link_data
    })
  end

  if step % 50 == 0 then
    local prog = step / total_steps
    local bar = terminal.progress_bar(prog, 16)
    io.write(string.format("\r[%s] Simulating: %d/%d steps (Recorded: %d frames)",
      bar, step, total_steps, #frames))
    io.flush()
  end
end

print("\n")

local metadata = {
  project = "fulcrum",
  scenario = "double_pendulum",
  integrator = "Forest-Ruth 4th-Order Symplectic",
  dt = dt * record_stride,
  num_frames = #frames,
  dof = model.dof,
  initial_energy = E0,
  created_at = os.date("%Y-%m-%d %H:%M:%S")
}

local output_path = "web/trajectory.json"
local json_str = json_exporter.export_trajectory(metadata, frames)
json_exporter.write_file(output_path, json_str)

local summary = {
  string.format("Simulation Time:       %.3f seconds", total_steps * dt),
  string.format("Recorded Frames:       %d frames at %.1f Hz", #frames, 1.0 / (dt * record_stride)),
  string.format("Integrator:            Forest-Ruth 4th-order symplectic"),
  string.format("Initial Total Energy:  %+.4f J", E0),
  string.format("Target File:           %s", output_path),
  string.format("File Size:             %.1f KB", #json_str / 1024.0)
}
print(terminal.draw_box("TRAJECTORY EXPORT TELEMETRY", summary, 76))
print(terminal.rgb_fg(0, 230, 118) .. terminal.BOLD .. "Web visualizer trajectory export completed successfully." .. terminal.RESET .. "\n")
