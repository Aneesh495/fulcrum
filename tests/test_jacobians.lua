local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local kinematics = require("fulcrum.dynamics.kinematics")
local jacobians = require("fulcrum.dynamics.jacobians")

runner.suite("Geometric Jacobians and IK")

runner.test("Jacobian matches forward velocity", function()
  local m = model.new()
  m:add_body({
    name = "link1",
    parent_id = 0,
    joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()),
    mass = 1.0,
    com = vec3.new(0.5, 0.0, 0.0)
  })
  m:add_body({
    name = "link2",
    parent_id = 1,
    joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.new(1.0, 0.0, 0.0)),
    mass = 1.0,
    com = vec3.new(0.5, 0.0, 0.0)
  })

  local q = { 0.5, -0.8 }
  local qd = { 1.2, -0.4 }

  -- Forward velocity of body 2
  local v, _, _ = kinematics.forward_velocities(m, q, qd)
  local v2_expected = v[2]

  -- Velocity from body-frame Jacobian
  local J2 = jacobians.geometric_jacobian(m, q, 2, false)
  local v2_jac = matrix.mul_vec(J2, qd)

  runner.assert_vector_near(v2_jac, v2_expected, 1e-10)
end)

runner.test("Nullspace projector property J * N = 0", function()
  local m = model.new()
  for i = 1, 3 do
    m:add_body({
      name = "link" .. i,
      parent_id = i - 1,
      joint = joint.revolute(vec3.unit_z()),
      X_T = spatial.transform(mat3.identity(), vec3.new(1.0, 0.0, 0.0)),
      mass = 1.0,
      com = vec3.new(0.5, 0.0, 0.0)
    })
  end

  local q = { 0.2, -0.5, 0.7 }
  local J = jacobians.geometric_jacobian(m, q, 3, true)
  local _, N = jacobians.nullspace_projector(J, 1e-8)

  local JN = J * N
  for i = 1, JN.rows do
    for j = 1, JN.cols do
      runner.assert_near(JN[i][j], 0.0, 1e-6)
    end
  end
end)
