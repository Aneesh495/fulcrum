local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local so3 = require("fulcrum.math.so3")
local se3 = require("fulcrum.math.se3")

runner.suite("SE(3) Lie Group Calculus")

runner.test("SE(3) multiplication and inverse", function()
  local R = so3.exp(vec3.new(0.5, -0.8, 1.2))
  local p = vec3.new(2.0, -3.5, 4.1)
  local T = se3.new(R, p)

  local T_inv = se3.inv(T)
  local prod = T * T_inv

  local I3 = mat3.identity()
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(prod.R[i][j], I3[i][j], 1e-10)
    end
  end
  runner.assert_near(prod.p.x, 0.0, 1e-10)
  runner.assert_near(prod.p.y, 0.0, 1e-10)
  runner.assert_near(prod.p.z, 0.0, 1e-10)
end)

runner.test("Exponential and logarithm roundtrip on SE(3)", function()
  local twist = { 1.5, -2.0, 3.2, 0.4, -0.6, 0.8 }
  local T = se3.exp(twist)
  local twist_rec = se3.log(T)

  for i = 1, 6 do
    runner.assert_near(twist_rec[i], twist[i], 1e-9)
  end
end)

runner.test("Point action transformation consistency", function()
  local T1 = se3.new(so3.exp(vec3.new(0.2, 0.3, -0.4)), vec3.new(1.0, 0.5, -2.0))
  local T2 = se3.new(so3.exp(vec3.new(-0.5, 0.1, 0.6)), vec3.new(-0.8, 1.2, 0.4))
  local pt = vec3.new(3.0, -1.0, 2.5)

  local pt_composed = se3.act(T1 * T2, pt)
  local pt_stepped = se3.act(T1, se3.act(T2, pt))

  runner.assert_near(pt_composed.x, pt_stepped.x, 1e-10)
  runner.assert_near(pt_composed.y, pt_stepped.y, 1e-10)
  runner.assert_near(pt_composed.z, pt_stepped.z, 1e-10)
end)

runner.test("Geodesic screw interpolation", function()
  local T1 = se3.identity()
  local twist = { 2.0, 0.0, 0.0, 0.0, 0.0, math.pi * 0.5 }
  local T2 = se3.exp(twist)

  local T_mid = se3.interpolate(T1, T2, 0.5)
  local twist_mid = se3.log(T_mid)

  runner.assert_near(twist_mid[1], 1.0, 1e-7)
  runner.assert_near(twist_mid[6], math.pi * 0.25, 1e-7)
end)
