local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local so3 = require("fulcrum.math.so3")

runner.suite("SO(3) Lie Group Calculus")

runner.test("Hat and vee isomorphism", function()
  local v = vec3.new(2.5, -1.8, 4.3)
  local V_hat = so3.hat(v)
  local v_rec = so3.vee(V_hat)

  runner.assert_near(v_rec.x, v.x, 1e-12)
  runner.assert_near(v_rec.y, v.y, 1e-12)
  runner.assert_near(v_rec.z, v.z, 1e-12)
end)

runner.test("Lie bracket and cross product equivalence", function()
  local u = vec3.new(1.1, 0.4, -2.2)
  local v = vec3.new(-0.8, 3.5, 1.0)

  local bracket = so3.bracket(so3.hat(u), so3.hat(v))
  local cross_hat = so3.hat(vec3.cross(u, v))

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(bracket[i][j], cross_hat[i][j], 1e-10)
    end
  end
end)

runner.test("Exponential and logarithm roundtrip", function()
  local phi = vec3.new(0.6, -1.2, 0.9)
  local R = so3.exp(phi)

  -- Check R in SO(3): R^T * R = I and det(R) = 1
  local RtR = mat3.transpose(R) * R
  local I = mat3.identity()
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(RtR[i][j], I[i][j], 1e-10)
    end
  end
  runner.assert_near(mat3.det(R), 1.0, 1e-10)

  local phi_rec = so3.log(R)
  runner.assert_near(phi_rec.x, phi.x, 1e-10)
  runner.assert_near(phi_rec.y, phi.y, 1e-10)
  runner.assert_near(phi_rec.z, phi.z, 1e-10)
end)

runner.test("Left and Right Jacobian inverses", function()
  local phi = vec3.new(0.4, 0.8, -0.5)
  local Jl = so3.left_jacobian(phi)
  local Jl_inv = so3.left_jacobian_inv(phi)

  local prod_l = Jl * Jl_inv
  local I = mat3.identity()
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(prod_l[i][j], I[i][j], 1e-10)
    end
  end

  local Jr = so3.right_jacobian(phi)
  local Jr_inv = so3.right_jacobian_inv(phi)
  local prod_r = Jr * Jr_inv
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(prod_r[i][j], I[i][j], 1e-10)
    end
  end
end)

runner.test("Adjoint action property", function()
  -- R * hat(phi) * R^T = hat(R * phi)
  local R = so3.exp(vec3.new(0.3, -0.7, 1.1))
  local phi = vec3.new(-1.5, 0.8, 2.2)

  local lhs = R * so3.hat(phi) * mat3.transpose(R)
  local rhs = so3.hat(R * phi)

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(lhs[i][j], rhs[i][j], 1e-10)
    end
  end
end)
