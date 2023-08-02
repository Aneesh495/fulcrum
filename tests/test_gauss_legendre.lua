local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local crba = require("fulcrum.dynamics.crba")
local gl = require("fulcrum.integrators.gauss_legendre")
local rkf45 = require("fulcrum.integrators.rkf45")

runner.suite("Gauss-Legendre and Adaptive RKF45")

runner.test("Gauss-Legendre 4th-order implicit integration", function()
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "link", parent_id = 0, joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(1.0, 0, 0),
    I_rot = mat3.from_diag(0.01, 0.01, 0.01)
  })

  local q = { 0.5 }
  local qd = { 0.0 }
  local dt = 0.01

  for step = 1, 100 do
    local acc
    q, qd, acc = gl.step4(m, q, qd, nil, dt)
  end

  local M = crba.inertia_matrix(m, q)
  local E = 0.5 * M[1][1] * qd[1] * qd[1] - 9.81 * math.sin(q[1])
  local E0 = -9.81 * math.sin(0.5)

  runner.assert_near(E, E0, 1e-4)
end)

runner.test("RKF45 adaptive step sizing", function()
  local m = model.new()
  m:add_body({
    name = "link", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(0.5, 0, 0)
  })

  local q = { 0.0 }
  local qd = { 1.0 }
  local dt = 0.05

  local q_next, v_next, dt_next, err = rkf45.step(m, q, qd, nil, dt, 1e-5)

  runner.assert_near(q_next[1], 0.05, 0.01)
  runner.assert_near(v_next[1], 1.0, 1e-5)
  runner.assert_true(dt_next > 0.0, "Next dt must be positive")
end)
