-- examples/03_newtons_cradle.lua
-- Demonstrates multi-rigid-body impact propagation across a 5-ball Newton's Cradle
-- using the Moreau-Jean velocity-impulse time-stepping scheme and Lemke LCP solver.

local vec3 = require("fulcrum.math.vec3")
local matrix = require("fulcrum.math.matrix")
local newtons_cradle = require("fulcrum.models.newtons_cradle")
local moreau_jean = require("fulcrum.contact.moreau_jean")
local terminal = require("fulcrum.viz.terminal")

print(terminal.CLEAR_SCREEN .. terminal.CURSOR_HOME)
print(terminal.rgb_fg(0, 229, 255) .. terminal.BOLD .. "=== FULCRUM: NEWTON'S CRADLE IMPACT PROPAGATION ===" .. terminal.RESET)
print(terminal.rgb_fg(139, 148, 158) .. "Moreau-Jean Time-Stepping with Lemke Complementarity Pivot Solver\n" .. terminal.RESET)

-- 5-ball cradle
local num_balls = 5
local model, meta = newtons_cradle.create({
  num_balls = num_balls,
  radius = 0.04,
  mass = 0.25,
  length = 0.45,
  gravity = vec3.new(0, 0, -9.81)
})

local L = meta.length
local R = meta.radius

-- Initial state: ball 1 displaced by -0.35 rad (pulled back to left)
local q = {}
local qd = {}
local tau = {}
for i = 1, num_balls do
  q[i] = (i == 1) and -0.35 or 0.0
  qd[i] = 0.0
  tau[i] = 0.0
end

local dt = 0.001
local total_steps = 1200
local print_interval = 20
local restitution = 0.99 -- Near-elastic impact transmission

local total_impacts = 0

local function render_ascii_cradle(q_curr)
  local width = 42
  local line = {}
  for i = 1, width do line[i] = " " end

  local center_col = math.floor(width / 2)
  for i = 1, num_balls do
    local nominal_x = (i - (num_balls + 1) * 0.5) * 4.0
    local disp_x = (q_curr[i] or 0.0) * 20.0
    local col = math.floor(center_col + nominal_x + disp_x + 0.5)
    if col >= 1 and col <= width then
      line[col] = (i == 1 or i == num_balls) and "O" or "o"
    end
  end

  return "[" .. table.concat(line) .. "]"
end

for step = 1, total_steps do
  local t = step * dt

  -- Detect contacts between adjacent balls (i and i+1)
  local active_contacts = {}
  for i = 1, num_balls - 1 do
    -- Gap: delta_x = 2*R + L*(sin(th_{i+1}) - sin(th_i)) - 2*R
    local th_a = q[i]
    local th_b = q[i + 1]
    local gap = L * (math.sin(th_b) - math.sin(th_a))

    if gap <= 1e-4 then
      -- Normal relative velocity: vn = L * (qd_{i+1}*cos(th_b) - qd_i*cos(th_a))
      local J_norm = {}
      for j = 1, num_balls do J_norm[j] = 0.0 end
      J_norm[i] = -L * math.cos(th_a)
      J_norm[i + 1] = L * math.cos(th_b)

      table.insert(active_contacts, {
        body_id = i,
        gap = gap,
        J_normal = J_norm
      })
    end
  end

  local P_n
  q, qd, P_n = moreau_jean.step(model, q, qd, tau, dt, active_contacts, restitution, 0.0, "lemke")

  if P_n and #P_n > 0 then
    local max_p = 0.0
    for _, p in ipairs(P_n) do
      if p > max_p then max_p = p end
    end
    if max_p > 1e-3 then
      total_impacts = total_impacts + 1
    end
  end

  if step % print_interval == 0 then
    local ascii_view = render_ascii_cradle(q)
    local prog = step / total_steps
    local bar = terminal.progress_bar(prog, 12)

    local v_str = {}
    for i = 1, num_balls do
      table.insert(v_str, string.format("%+.2f", qd[i]))
    end

    io.write(string.format("\r[%s] t=%.3fs %s | v=[%s] | Impacts: %d",
      bar, t, ascii_view, table.concat(v_str, " "), total_impacts))
    io.flush()
  end
end

print("\n")

local summary = {
  string.format("Simulation Time:    %.3f seconds (%d steps at dt = %.4f s)", total_steps * dt, total_steps, dt),
  string.format("Number of Spheres:  %d identical steel balls", num_balls),
  string.format("Restitution (e):    %.2f (near-elastic)", restitution),
  string.format("Impact Events:      %d resolved via Moreau-Jean Lemke LCP", total_impacts),
  string.format("Final Ball 1 Speed: %+.4f rad/s", qd[1]),
  string.format("Final Ball 5 Speed: %+.4f rad/s (wave propagated to opposite end)", qd[num_balls])
}
print(terminal.draw_box("IMPACT PROPAGATION TELEMETRY", summary, 74))
print(terminal.rgb_fg(0, 230, 118) .. terminal.BOLD .. "Newton's Cradle contact wave propagation completed successfully." .. terminal.RESET .. "\n")
