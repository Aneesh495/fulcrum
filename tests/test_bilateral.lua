local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local bilateral = require("fulcrum.contact.bilateral")
local moreau_jean = require("fulcrum.contact.moreau_jean")

runner.suite("Bilateral Constraints and Moreau-Jean Time-Stepping")

runner.test("Bilateral constraint acceleration and Lagrange multipliers", function()
  -- Cart on rail with constraint x = 0 (prismatic joint along X, constraint forces q[1] = 0)
  local m = model.new()
  m:add_body({
    name = "cart", parent_id = 0, joint = joint.prismatic(vec3.unit_x()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 5.0, com = vec3.zero()
  })

  -- Constraint: g(q) = q[1] = 0 -> G = [ 1 ]
  local G = matrix.from_table({ { 1.0 } })
  local g_val = { 0.0 }
  local g_dot_val = { 0.0 }
  local tau = { 50.0 } -- apply 50 N force along X

  local qdd, lambda = bilateral.solve_constrained_acceleration(m, { 0.0 }, { 0.0 }, tau, G, g_val, g_dot_val, 10.0, 10.0)

  -- Cart acceleration must be zero because constraint resists it
  runner.assert_near(qdd[1], 0.0, 1e-6)
  -- Constraint force lambda must balance applied force 50 N
  runner.assert_near(lambda[1], 50.0, 1e-4)
end)

runner.test("Moreau-Jean contact impact resolution", function()
  -- Falling mass hitting ground
  local m = model.new({ gravity = vec3.new(0, 0, -9.81) })
  m:add_body({
    name = "ball", parent_id = 0, joint = joint.prismatic(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 2.0, com = vec3.zero()
  })

  -- Downward velocity: qd = -5.0 m/s
  local q = { 0.0 }
  local qd = { -5.0 }
  local dt = 0.01

  -- Contact constraint: normal along +Z -> J_normal = [ 1 ]
  local contacts = {
    {
      body_id = 1,
      normal = vec3.unit_z(),
      J_normal = { 1.0 }
    }
  }

  -- Inelastic impact: restitution e = 0
  local q_next, v_next, P_n = moreau_jean.step(m, q, qd, nil, dt, contacts, 0.0)

  -- Post-impact velocity should be non-negative (stopped at surface)
  runner.assert_near(v_next[1], 0.0, 1e-4)
  -- Contact impulse must be positive to cancel momentum (m * delta_v = 2.0 * 5.0 = 10.0 Ns approx)
  runner.assert_true(P_n[1] > 9.0, "Impulse must counteract initial momentum")

  -- Elastic impact: restitution e = 0.8
  local _, v_bounce, _ = moreau_jean.step(m, q, qd, nil, dt, contacts, 0.8)
  -- Rebound velocity: v_next = -e * qd = 0.8 * 5.0 = 4.0 m/s approx
  runner.assert_near(v_bounce[1], 4.0, 0.1)
end)
