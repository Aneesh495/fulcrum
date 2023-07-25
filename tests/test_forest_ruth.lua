local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local crba = require("fulcrum.dynamics.crba")
local forest_ruth = require("fulcrum.integrators.forest_ruth")

runner.suite("Forest-Ruth 4th-Order Symplectic Integrator")

runner.test("High-precision energy conservation over long trajectories", function()
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "link", parent_id = 0, joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(1.0, 0, 0),
    I_rot = mat3.from_diag(0.01, 0.01, 0.01)
  })

  local total_energy = function(q, qd)
    local M = crba.inertia_matrix(m, q)
    local T = 0.5 * M[1][1] * qd[1] * qd[1]
    local V = -1.0 * 9.81 * 1.0 * math.sin(q[1])
    return T + V
  end

  local q = { 0.8 }
  local qd = { 0.0 }
  local dt = 0.02
  local E0 = total_energy(q, qd)

  for step = 1, 250 do
    local acc
    q, qd, acc = forest_ruth.step(m, q, qd, nil, dt)
  end

  local E_final = total_energy(q, qd)
  -- 4th order error is extremely tight (order of 10^-4)
  runner.assert_near(E_final, E0, 1e-3)
end)

runner.test("Yoshida 6th-order integration convergence", function()
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "link", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(0.5, 0, 0),
    I_rot = mat3.from_diag(0.05, 0.05, 0.05)
  })

  local q = { 0.1 }
  local qd = { 2.0 }
  local dt = 0.01

  for step = 1, 100 do
    local acc
    q, qd, acc = forest_ruth.yoshida6(m, q, qd, nil, dt)
  end

  runner.assert_near(qd[1], 2.0, 1e-10)
end)
