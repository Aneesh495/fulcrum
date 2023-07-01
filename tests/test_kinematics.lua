local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local kinematics = require("fulcrum.dynamics.kinematics")

runner.suite("Kinematics and Tree Propagation")

runner.test("Two-link planar arm forward kinematics", function()
  local m = model.new()

  -- Link 1: 1m length along X, rotates around Z
  local l1 = 1.0
  m:add_body({
    name = "link1",
    parent_id = 0,
    joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()),
    mass = 2.0,
    com = vec3.new(l1 * 0.5, 0.0, 0.0)
  })

  -- Link 2: 1m length along X, connected at end of link 1
  local l2 = 1.0
  m:add_body({
    name = "link2",
    parent_id = 1,
    joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.new(l1, 0.0, 0.0)),
    mass = 1.5,
    com = vec3.new(l2 * 0.5, 0.0, 0.0)
  })

  runner.assert_equal(m.dof, 2)

  -- Position at q = [0, 0]
  local q0 = { 0.0, 0.0 }
  local _, X_world = kinematics.forward_kinematics(m, q0)
  runner.assert_equal(#X_world, 2)

  -- Forward velocities at q = [0, 0], q_dot = [1.0, 2.0]
  local qd = { 1.0, 2.0 }
  local v = kinematics.forward_velocities(m, q0, qd)
  runner.assert_equal(#v, 2)

  -- Link 1 angular velocity around Z is 1.0
  runner.assert_near(v[1][3], 1.0, 1e-10)
  -- Link 2 angular velocity around Z is 1.0 + 2.0 = 3.0
  runner.assert_near(v[2][3], 3.0, 1e-10)
end)

runner.test("Center of mass computation", function()
  local m = model.new()
  m:add_body({
    name = "body1",
    parent_id = 0,
    joint = joint.prismatic(vec3.unit_x()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()),
    mass = 10.0,
    com = vec3.new(0.0, 0.0, 0.0)
  })

  local q = { 5.0 } -- translated by 5 along X
  local com, total_m = kinematics.center_of_mass(m, q)

  runner.assert_near(total_m, 10.0, 1e-10)
  runner.assert_near(com.x, 5.0, 1e-7)
  runner.assert_near(com.y, 0.0, 1e-7)
  runner.assert_near(com.z, 0.0, 1e-7)
end)
