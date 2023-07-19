local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local jacobians = require("fulcrum.dynamics.jacobians")
local crba = require("fulcrum.dynamics.crba")
local op_space = require("fulcrum.dynamics.operational_space")

runner.suite("Operational Space Formulation")

runner.test("Dynamically consistent inverse and decoupling properties", function()
  -- 3-link planar manipulator (redundant 3 DOF for 2D position task)
  local m = model.new()
  for i = 1, 3 do
    m:add_body({
      name = "link" .. i,
      parent_id = i - 1,
      joint = joint.revolute(vec3.unit_z()),
      X_T = spatial.transform(mat3.identity(), vec3.new(1.0, 0.0, 0.0)),
      mass = 1.5,
      com = vec3.new(0.5, 0.0, 0.0),
      I_rot = mat3.from_diag(0.1, 0.1, 0.1)
    })
  end

  local q = { 0.4, -0.8, 0.5 }
  local J_full = jacobians.geometric_jacobian(m, q, 3, true)
  -- 2D planar position task (rows 4 and 5 in world frame)
  local J_task = matrix.new(2, 3)
  for c = 1, 3 do
    J_task[1][c] = J_full[4][c]
    J_task[2][c] = J_full[5][c]
  end

  local M = crba.inertia_matrix(m, q)
  local Lambda, M_inv = op_space.inertia(M, J_task)
  local J_sharp = op_space.dynamically_consistent_inverse(M_inv, J_task, Lambda)

  -- Property 1: J * J^# = I_task
  local JJ_sharp = J_task * J_sharp
  local I2 = matrix.eye(2)
  for i = 1, 2 do
    for j = 1, 2 do
      runner.assert_near(JJ_sharp[i][j], I2[i][j], 1e-8)
    end
  end

  -- Property 2: Nullspace projection orthogonality in task space
  local N_T = op_space.nullspace_projector(J_task, J_sharp)
  local tau_null = { 5.0, -3.0, 2.0 }
  local tau_proj = matrix.mul_vec(N_T, tau_null)

  -- Effective task force produced by tau_proj: F_eff = (J^#)^T * tau_proj == 0
  local J_sharp_t = matrix.transpose(J_sharp)
  local F_eff = matrix.mul_vec(J_sharp_t, tau_proj)
  runner.assert_near(F_eff[1], 0.0, 1e-8)
  runner.assert_near(F_eff[2], 0.0, 1e-8)
end)
