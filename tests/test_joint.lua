local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local matrix = require("fulcrum.math.matrix")
local joint = require("fulcrum.dynamics.joint")

runner.suite("Articulated Joint Models")

runner.test("Revolute joint subspace and transform", function()
  local j = joint.revolute(vec3.unit_z(), "elbow")
  runner.assert_equal(j.dof, 1)

  local S = j.subspace(0.0)
  runner.assert_equal(S.rows, 6)
  runner.assert_equal(S.cols, 1)
  runner.assert_near(S[1][1], 0.0)
  runner.assert_near(S[2][1], 0.0)
  runner.assert_near(S[3][1], 1.0)
  runner.assert_near(S[4][1], 0.0)

  local X_J = j.transform(math.pi * 0.5)
  -- 90 deg rotation around Z
  runner.assert_near(X_J[1][2], -1.0, 1e-10)
  runner.assert_near(X_J[2][1], 1.0, 1e-10)
end)

runner.test("Prismatic joint subspace and translation", function()
  local j = joint.prismatic(vec3.unit_x(), "slider")
  runner.assert_equal(j.dof, 1)

  local S = j.subspace(0.0)
  runner.assert_near(S[4][1], 1.0)
  runner.assert_near(S[5][1], 0.0)
  runner.assert_near(S[6][1], 0.0)

  local X_J = j.transform(2.5)
  -- Bottom-left 3x3 of spatial transform has -R*skew(p)
  -- p = [2.5, 0, 0] -> skew(p) = [0, 0, 0; 0, 0, -2.5; 0, 2.5, 0]
  -- -R*skew(p) = [0, 0, 0; 0, 0, 2.5; 0, -2.5, 0]
  runner.assert_near(X_J[5][3], 2.5, 1e-10)
  runner.assert_near(X_J[6][2], -2.5, 1e-10)
end)

runner.test("Helical joint coupled motion", function()
  local pitch = 0.05
  local j = joint.helical(vec3.unit_z(), pitch, "screw")
  local S = j.subspace(0.0)

  runner.assert_near(S[3][1], 1.0)
  runner.assert_near(S[6][1], pitch)
end)
