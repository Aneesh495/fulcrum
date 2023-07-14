local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")
local spatial = require("fulcrum.math.spatial")
local joint = require("fulcrum.dynamics.joint")
local model = require("fulcrum.dynamics.model")
local rnea = require("fulcrum.dynamics.rnea")
local crba = require("fulcrum.dynamics.crba")

runner.suite("Composite Rigid Body Dynamics (CRBA)")

runner.test("Symmetry and positive-definiteness of M(q)", function()
  local m = model.new()
  m:add_body({
    name = "link1", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 2.0, com = vec3.new(0.5, 0, 0),
    I_rot = mat3.from_diag(0.2, 0.2, 0.2)
  })
  m:add_body({
    name = "link2", parent_id = 1, joint = joint.revolute(vec3.unit_y()),
    X_T = spatial.transform(mat3.identity(), vec3.new(1.0, 0, 0)), mass = 1.5, com = vec3.new(0.4, 0, 0),
    I_rot = mat3.from_diag(0.1, 0.1, 0.1)
  })
  m:add_body({
    name = "link3", parent_id = 2, joint = joint.revolute(vec3.unit_x()),
    X_T = spatial.transform(mat3.identity(), vec3.new(0.8, 0, 0)), mass = 1.0, com = vec3.new(0.3, 0, 0),
    I_rot = mat3.from_diag(0.05, 0.05, 0.05)
  })

  local q = { 0.4, -0.7, 1.1 }
  local M = crba.inertia_matrix(m, q)

  -- Check exact symmetry
  for i = 1, m.dof do
    for j = 1, m.dof do
      runner.assert_near(M[i][j], M[j][i], 1e-12)
    end
  end

  -- Check positive-definiteness via Cholesky factorisation
  local L = decomp.cholesky(M)
  runner.assert_equal(L.rows, m.dof)
  for i = 1, m.dof do
    runner.assert_true(L[i][i] > 1e-4, "Diagonal pivots must be positive")
  end
end)

runner.test("CRBA matches RNEA M(q) * q_ddot identity", function()
  -- The fundamental multibody identity:
  -- M(q) * q_ddot == RNEA(q, 0, q_ddot) - RNEA(q, 0, 0)
  local m = model.new()
  m:add_body({
    name = "link1", parent_id = 0, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.zero()), mass = 3.0, com = vec3.new(0.6, 0.1, 0),
    I_rot = mat3.from_diag(0.3, 0.3, 0.3)
  })
  m:add_body({
    name = "link2", parent_id = 1, joint = joint.revolute(vec3.unit_z()),
    X_T = spatial.transform(mat3.identity(), vec3.new(1.2, 0, 0)), mass = 2.0, com = vec3.new(0.5, 0, 0),
    I_rot = mat3.from_diag(0.2, 0.2, 0.2)
  })

  local q = { -0.5, 0.9 }
  local qd_zero = { 0.0, 0.0 }
  local qdd = { 1.8, -2.4 }

  local M = crba.inertia_matrix(m, q)
  local M_qdd = matrix.mul_vec(M, qdd)

  local tau_acc = rnea.inverse_dynamics(m, q, qd_zero, qdd)
  local tau_zero = rnea.inverse_dynamics(m, q, qd_zero, qd_zero)

  local rnea_M_qdd = {
    tau_acc[1] - tau_zero[1],
    tau_acc[2] - tau_zero[2]
  }

  runner.assert_near(M_qdd[1], rnea_M_qdd[1], 1e-9)
  runner.assert_near(M_qdd[2], rnea_M_qdd[2], 1e-9)
end)
