local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local rnea = require("fulcrum.dynamics.rnea")

runner.suite("Recursive Newton-Euler Dynamics (RNEA)")

runner.test("Single pendulum static gravity torque", function()
  -- Gravity pointing downwards along -Z: [0, 0, -9.81]
  local g = 9.81
  local m_mass = 2.0
  local l = 1.0

  local m = model.new({ gravity = vec3.new(0, 0, -g) })
  -- Rotate around Y axis, link extends along X axis
  m:add_body({
    name = "link",
    parent_id = 0,
    joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()),
    mass = m_mass,
    com = vec3.new(l, 0.0, 0.0),
    I_rot = mat3.from_diag(0.1, 0.1, 0.1)
  })

  -- When q = 0 (horizontal), gravity produces torque tau_y = -m * g * l (or +m*g*l)
  local q0 = { 0.0 }
  local tau_g0 = rnea.gravity_terms(m, q0)
  runner.assert_near(math.abs(tau_g0[1]), m_mass * g * l, 1e-10)

  -- When q = 90 deg (hanging downwards), torque should be 0
  local q_down = { math.pi * 0.5 }
  local tau_down = rnea.gravity_terms(m, q_down)
  runner.assert_near(tau_down[1], 0.0, 1e-10)
end)

runner.test("Linearity with respect to acceleration", function()
  local m = model.new()
  m:add_body({
    name = "link1", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 1.0, com = vec3.new(0.5, 0, 0)
  })
  m:add_body({
    name = "link2", parent_id = 1, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.new(1.0, 0, 0)), mass = 1.0, com = vec3.new(0.5, 0, 0)
  })

  local q = { 0.3, -0.6 }
  local qd = { 0.0, 0.0 }

  local qdd1 = { 1.0, 0.0 }
  local qdd2 = { 0.0, 1.0 }
  local qdd_comb = { 2.5, -1.8 }

  local tau0 = rnea.inverse_dynamics(m, q, qd, { 0.0, 0.0 })
  local tau1 = rnea.inverse_dynamics(m, q, qd, qdd1)
  local tau2 = rnea.inverse_dynamics(m, q, qd, qdd2)
  local tau_comb = rnea.inverse_dynamics(m, q, qd, qdd_comb)

  -- (tau1 - tau0)*2.5 + (tau2 - tau0)*(-1.8) + tau0 must equal tau_comb
  local exp1 = (tau1[1] - tau0[1]) * 2.5 + (tau2[1] - tau0[1]) * (-1.8) + tau0[1]
  local exp2 = (tau1[2] - tau0[2]) * 2.5 + (tau2[2] - tau0[2]) * (-1.8) + tau0[2]

  runner.assert_near(tau_comb[1], exp1, 1e-10)
  runner.assert_near(tau_comb[2], exp2, 1e-10)
end)
