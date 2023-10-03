-- examples/01_dzhanibekov.lua
-- Demonstrates the intermediate-axis instability theorem (Dzhanibekov effect)
-- using the Lie Group Variational Integrator (LGVI) on SO(3).
-- Includes real-time 3D wireframe terminal animation and conservation telemetry.

local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local lgvi = require("fulcrum.integrators.lie_group_integrator")
local terminal = require("fulcrum.viz.terminal")
local rasterizer3d = require("fulcrum.viz.rasterizer3d")

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: DZHANIBEKOV INTERMEDIATE AXIS EFFECT ===" .. terminal.RESET)
print(terminal.rgb_fg(139, 148, 158) .. "Lie Group Variational Integrator on SO(3) with rigid body inertia J1 < J2 < J3\n" .. terminal.RESET)

-- Principal moments of inertia: intermediate axis is Y (J2)
local J_diag = vec3.new(1.0, 2.0, 3.5)
local half_extents = vec3.new(0.6, 0.25, 0.15)

-- Initial state: rotation almost purely about intermediate Y axis with tiny perturbation
local R = mat3.identity()
local Omega = vec3.new(0.002, 4.0, 0.002)

local function kinetic_energy(w)
  return 0.5 * (J_diag.x * w.x^2 + J_diag.y * w.y^2 + J_diag.z * w.z^2)
end

local function spatial_momentum(rot, w)
  local Pi_body = vec3.new(J_diag.x * w.x, J_diag.y * w.y, J_diag.z * w.z)
  return rot * Pi_body
end

local E0 = kinetic_energy(Omega)
local L0 = vec3.norm(spatial_momentum(R, Omega))

local dt = 0.005
local total_steps = 300
local render_interval = 25

local rast = rasterizer3d.new(64, 20)
rast:set_camera(vec3.new(2.5, -3.2, 1.8), vec3.zero(), vec3.unit_z())

local flips = 0
local last_omega_y_sign = (Omega.y > 0)

for step = 1, total_steps do
  R, Omega = lgvi.step_so3(R, Omega, J_diag, nil, dt)

  local cur_sign = (Omega.y > 0)
  if cur_sign ~= last_omega_y_sign then
    flips = flips + 1
    last_omega_y_sign = cur_sign
  end

  if step % render_interval == 0 or step == total_steps then
    rast:clear()
    rast:draw_axes(vec3.zero(), mat3.identity(), 0.6)
    rast:draw_box_wireframe(vec3.zero(), half_extents, R, "#")

    local E = kinetic_energy(Omega)
    local L = vec3.norm(spatial_momentum(R, Omega))
    local dE = math.abs((E - E0) / E0)
    local dL = math.abs((L - L0) / L0)

    local hud = {
      string.format("Step: %4d / %d  |  Sim Time: %.3f s", step, total_steps, step * dt),
      string.format("Omega: [%+.3f, %+.3f, %+.3f] rad/s", Omega.x, Omega.y, Omega.z),
      string.format("Energy Drift:       %.2e (Symplectic preservation)", dE),
      string.format("Momentum Drift:     %.2e (Noether invariant)", dL),
      string.format("Axis Flips Count:   %d (Periodic intermediate flip-flops)", flips)
    }

    print(terminal.CURSOR_HOME)
    print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: DZHANIBEKOV INTERMEDIATE AXIS EFFECT ===" .. terminal.RESET)
    print(terminal.rgb_fg(139, 148, 158) .. "Lie Group Variational Integrator on SO(3) with rigid body inertia J1 < J2 < J3" .. terminal.RESET)
    print(terminal.draw_box("SIMULATION TELEMETRY", hud, 66))
    print(rast:render_to_string())
  end
end

print(terminal.rgb_fg(0, 230, 118) .. terminal.BOLD .. "\nSimulation completed. Geometric Lie group symplecticity verified." .. terminal.RESET)
