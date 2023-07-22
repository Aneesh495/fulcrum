local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local crba = require("fulcrum.dynamics.crba")
local rnea = require("fulcrum.dynamics.rnea")
local symplectic = require("fulcrum.integrators.symplectic")

runner.suite("Symplectic Integrators")

runner.test("Energy conservation of nonlinear pendulum over 500 steps", function()
  -- Simple undamped pendulum in zero gravity with spring potential or standard gravity
  local m_mass = 1.0
  local l = 1.0
  local g = 9.81

  local m = model.new({ gravity = vec3.new(0, 0, -g) })
  m:add_body({
    name = "bob", parent_id = 0, joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = m_mass, com = vec3.new(l, 0, 0),
    I_rot = mat3.from_diag(0.01, 0.01, 0.01)
  })

  local total_energy = function(q, qd)
    local M = crba.inertia_matrix(m, q)
    local T = 0.5 * M[1][1] * qd[1] * qd[1]
    -- Potential energy: V = -m * g * l * sin(q[1])
    local V = -m_mass * g * l * math.sin(q[1])
    return T + V
  end

  local q = { 0.5 }
  local qd = { 0.0 }
  local dt = 0.01
  local E0 = total_energy(q, qd)

  for step = 1, 500 do
    local acc
    q, qd, acc = symplectic.verlet(m, q, qd, nil, dt)
  end

  local E_final = total_energy(q, qd)
  -- Symplectic Verlet has bounded energy oscillation without secular drift
  runner.assert_near(E_final, E0, 0.05)
end)

runner.test("Implicit midpoint quadratic invariant preservation", function()
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "link", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(0.5, 0, 0),
    I_rot = mat3.from_diag(0.05, 0.05, 0.05)
  })

  local q = { 0.2 }
  local qd = { 1.5 }
  local dt = 0.005

  for step = 1, 100 do
    local acc
    q, qd, acc = symplectic.implicit_midpoint(m, q, qd, nil, dt)
  end

  runner.assert_near(qd[1], 1.5, 1e-10)
end)
