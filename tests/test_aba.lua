local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local rnea = require("fulcrum.dynamics.rnea")
local aba = require("fulcrum.dynamics.aba")

runner.suite("Articulated Body Forward Dynamics (ABA)")

runner.test("ABA forward dynamics inverted by RNEA: RNEA(q, qd, ABA(q, qd, tau)) == tau", function()
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "link1", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 2.5, com = vec3.new(0.4, 0, 0),
    I_rot = mat3.from_diag(0.15, 0.15, 0.15)
  })
  m:add_body({
    name = "link2", parent_id = 1, joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.new(0.8, 0, 0)), mass = 1.8, com = vec3.new(0.3, 0, 0),
    I_rot = mat3.from_diag(0.1, 0.1, 0.1)
  })
  m:add_body({
    name = "link3", parent_id = 2, joint = joint.revolute(vec3.unit_x()),
    X_T = spatial.transform(mat3.identity(), vec3.new(0.6, 0, 0)), mass = 1.2, com = vec3.new(0.2, 0, 0),
    I_rot = mat3.from_diag(0.05, 0.05, 0.05)
  })

  local q = { 0.7, -1.1, 0.4 }
  local qd = { -1.5, 0.8, 2.1 }
  local tau = { 15.0, -8.5, 4.2 }

  -- Forward dynamics: compute joint accelerations
  local qdd = aba.forward_dynamics(m, q, qd, tau)
  runner.assert_equal(#qdd, 3)

  -- Inverse dynamics: recover applied torques from computed accelerations
  local tau_recovered = rnea.inverse_dynamics(m, q, qd, qdd)

  for i = 1, 3 do
    runner.assert_near(tau_recovered[i], tau[i], 1e-9)
  end
end)
