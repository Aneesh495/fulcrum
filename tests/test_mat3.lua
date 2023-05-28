local runner = require("tests.test_framework")
local mat3 = require("fulcrum.math.mat3")
local vec3 = require("fulcrum.math.vec3")

runner.suite("3x3 Matrix Algebra")

runner.test("Identity and basic properties", function()
  local I = mat3.identity()
  runner.assert_near(mat3.det(I), 1.0)
  runner.assert_near(mat3.trace(I), 3.0)

  local v = vec3.new(3.5, -2.1, 8.4)
  local Iv = I * v
  runner.assert_near(Iv.x, v.x)
  runner.assert_near(Iv.y, v.y)
  runner.assert_near(Iv.z, v.z)
end)

runner.test("Matrix multiplication and associativity", function()
  local A = mat3.new(
    1, 2, 3,
    0, 1, 4,
    5, 6, 0
  )
  local B = mat3.new(
    2, 0, -1,
    1, 3, 2,
    0, -2, 1
  )
  local C = mat3.new(
    1, 1, 0,
    0, 2, 1,
    1, 0, 1
  )

  local AB_C = (A * B) * C
  local A_BC = A * (B * C)

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(AB_C[i][j], A_BC[i][j], 1e-10)
    end
  end
end)

runner.test("Matrix inverse and determinant", function()
  local A = mat3.new(
    2, -1, 0,
    -1, 2, -1,
    0, -1, 2
  )
  local detA = mat3.det(A)
  runner.assert_near(detA, 4.0)

  local Ainv = mat3.inv(A)
  local prod = A * Ainv
  local I = mat3.identity()

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(prod[i][j], I[i][j], 1e-10)
    end
  end
end)

runner.test("Axis-angle matrix generation", function()
  local axis = vec3.unit(vec3.new(1.0, 1.0, 1.0))
  local angle = 2.0 * math.pi / 3.0 -- 120 degrees

  local R = mat3.from_axis_angle(axis, angle)
  -- 120 deg around [1,1,1] permutes x -> y -> z -> x
  local vx = vec3.unit_x()
  local vy = vec3.unit_y()
  local vz = vec3.unit_z()

  local Rvx = R * vx
  runner.assert_near(Rvx.x, vy.x, 1e-10)
  runner.assert_near(Rvx.y, vy.y, 1e-10)
  runner.assert_near(Rvx.z, vy.z, 1e-10)

  local rec_axis, rec_angle = mat3.to_axis_angle(R)
  runner.assert_near(rec_angle, angle, 1e-7)
  runner.assert_near(rec_axis.x, axis.x, 1e-7)
  runner.assert_near(rec_axis.y, axis.y, 1e-7)
  runner.assert_near(rec_axis.z, axis.z, 1e-7)
end)

runner.test("Gram-Schmidt orthonormalization", function()
  -- Distorted non-orthogonal frame
  local M = mat3.new(
    1.05, 0.1, -0.05,
    -0.08, 0.98, 0.12,
    0.04, -0.09, 1.03
  )
  local R = mat3.orthonormalize(M)

  -- Check R^T * R = I
  local RtR = mat3.transpose(R) * R
  local I = mat3.identity()

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(RtR[i][j], I[i][j], 1e-10)
    end
  end
  runner.assert_near(mat3.det(R), 1.0, 1e-10)
end)
