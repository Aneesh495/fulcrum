local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local models = require("fulcrum.models")
local crba = require("fulcrum.dynamics.crba")
local rnea = require("fulcrum.dynamics.rnea")
local aba = require("fulcrum.dynamics.aba")
local jacobians = require("fulcrum.dynamics.jacobians")

runner.suite("Benchmark Physics and Articulated Robot Models")

runner.test("Cartpole 2-DOF underactuated system", function()
  local model, meta = models.cartpole.create()
  runner.assert_equal(model.num_bodies, 2, "Cartpole has 2 bodies")
  runner.assert_equal(model.dof, 2, "Cartpole has 2 DOF")

  local q = { 0.5, 0.2 } -- cart at 0.5m, pole at 0.2 rad
  local M = crba.mass_matrix(model, q)
  runner.assert_equal(M.rows, 2, "Mass matrix is 2x2")
  runner.assert_true(M[1][1] > 0.0, "Cart effective mass > 0")
  runner.assert_true(M[2][2] > 0.0, "Pole effective inertia > 0")
  runner.assert_near(M[1][2], M[2][1], 1e-6, "Mass matrix is symmetric")
end)

runner.test("Double and triple chaotic pendulums", function()
  local dp, _ = models.pendulums.double_pendulum()
  runner.assert_equal(dp.dof, 2, "Double pendulum has 2 DOF")

  local q2 = { 0.3, -0.4 }
  local qd2 = { 0.1, 0.5 }
  local tau2 = { 0.0, 0.0 }
  local qdd2 = aba.forward_dynamics(dp, q2, qd2, tau2)
  local tau_check2 = rnea.inverse_dynamics(dp, q2, qd2, qdd2)
  runner.assert_near(tau_check2[1], 0.0, 1e-6, "Double pendulum ABA/RNEA consistency 1")
  runner.assert_near(tau_check2[2], 0.0, 1e-6, "Double pendulum ABA/RNEA consistency 2")

  local tp, _ = models.pendulums.triple_pendulum()
  runner.assert_equal(tp.dof, 3, "Triple pendulum has 3 DOF")
end)

runner.test("Dzhanibekov intermediate axis rigid body", function()
  local model, meta = models.dzhanibekov.create({ I1 = 1.0, I2 = 2.0, I3 = 3.0 })
  runner.assert_equal(model.num_bodies, 1, "Dzhanibekov is a single rigid body")
  runner.assert_equal(model.dof, 6, "Floating base has 6 DOF")
  runner.assert_true(meta.I1 < meta.I2 and meta.I2 < meta.I3, "Intermediate axis inequality I1 < I2 < I3")
end)

runner.test("Newton's Cradle 5-sphere impact array", function()
  local model, meta = models.newtons_cradle.create({ num_balls = 5 })
  runner.assert_equal(model.num_bodies, 5, "5 suspended pendulums")
  runner.assert_equal(#meta.shapes, 5, "5 collision sphere shapes created")
end)

runner.test("7-DOF Anthropomorphic Manipulator and Kinematic Jacobian", function()
  local model, meta = models.manipulator_7dof.create()
  runner.assert_equal(model.dof, 7, "Manipulator has 7 DOF")

  local q = { 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, -0.7 }
  local J_ee = jacobians.body_jacobian(model, q, meta.end_effector_id)
  runner.assert_equal(J_ee.rows, 6, "Jacobian has 6 spatial rows")
  runner.assert_equal(J_ee.cols, 7, "Jacobian has 7 joint columns")
end)

runner.test("Biped walker and Quadruped 12-DOF mobile platforms", function()
  local biped_fixed, b_meta = models.biped.create({ floating = false })
  runner.assert_equal(biped_fixed.dof, 6, "Fixed-base biped has 6 leg joints")
  runner.assert_true(b_meta.left_foot_id ~= nil, "Left foot ID mapped")
  runner.assert_true(b_meta.right_foot_id ~= nil, "Right foot ID mapped")

  local quad_fixed, q_meta = models.quadruped.create({ floating = false })
  runner.assert_equal(quad_fixed.dof, 12, "Fixed-base quadruped has 12 leg joints")
  runner.assert_true(q_meta.foot_ids.fl ~= nil, "Front-left foot mapped")
  runner.assert_true(q_meta.foot_ids.rr ~= nil, "Rear-right foot mapped")
end)
