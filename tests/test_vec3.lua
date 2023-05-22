local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")

runner.suite("Vector 3D Algebra")

runner.test("Vector creation and components", function()
  local v = vec3.new(1.0, 2.0, 3.0)
  runner.assert_equal(v.x, 1.0)
  runner.assert_equal(v.y, 2.0)
  runner.assert_equal(v.z, 3.0)
end)

runner.test("Vector arithmetic operators", function()
  local a = vec3.new(1.0, 2.0, 3.0)
  local b = vec3.new(4.0, 5.0, 6.0)

  local sum = a + b
  runner.assert_near(sum.x, 5.0)
  runner.assert_near(sum.y, 7.0)
  runner.assert_near(sum.z, 9.0)

  local diff = b - a
  runner.assert_near(diff.x, 3.0)
  runner.assert_near(diff.y, 3.0)
  runner.assert_near(diff.z, 3.0)

  local scaled = a * 2.5
  runner.assert_near(scaled.x, 2.5)
  runner.assert_near(scaled.y, 5.0)
  runner.assert_near(scaled.z, 7.5)

  local div = scaled / 2.5
  runner.assert_near(div.x, a.x)
  runner.assert_near(div.y, a.y)
  runner.assert_near(div.z, a.z)

  local neg = -a
  runner.assert_near(neg.x, -1.0)
  runner.assert_near(neg.y, -2.0)
  runner.assert_near(neg.z, -3.0)
end)

runner.test("Dot and cross product", function()
  local i = vec3.unit_x()
  local j = vec3.unit_y()
  local k = vec3.unit_z()

  runner.assert_near(vec3.dot(i, j), 0.0)
  runner.assert_near(vec3.dot(i, i), 1.0)

  local cross_ij = vec3.cross(i, j)
  runner.assert_near(cross_ij.x, k.x)
  runner.assert_near(cross_ij.y, k.y)
  runner.assert_near(cross_ij.z, k.z)

  local cross_jk = vec3.cross(j, k)
  runner.assert_near(cross_jk.x, i.x)
  runner.assert_near(cross_jk.y, i.y)
  runner.assert_near(cross_jk.z, i.z)

  -- Anti-symmetry: a x b = -(b x a)
  local a = vec3.new(1.2, -3.4, 5.6)
  local b = vec3.new(-7.8, 9.0, 2.1)
  local c1 = vec3.cross(a, b)
  local c2 = -vec3.cross(b, a)
  runner.assert_near(c1.x, c2.x)
  runner.assert_near(c1.y, c2.y)
  runner.assert_near(c1.z, c2.z)

  -- Lagrange identity: |a x b|^2 = |a|^2 |b|^2 - (a . b)^2
  local lhs = vec3.norm_sq(c1)
  local rhs = vec3.norm_sq(a) * vec3.norm_sq(b) - (vec3.dot(a, b))^2
  runner.assert_near(lhs, rhs, 1e-10)
end)

runner.test("Projection, rejection, and reflection", function()
  local v = vec3.new(3.0, 4.0, 5.0)
  local n = vec3.unit_z()

  local p = vec3.project(v, n)
  runner.assert_near(p.x, 0.0)
  runner.assert_near(p.y, 0.0)
  runner.assert_near(p.z, 5.0)

  local r = vec3.reject(v, n)
  runner.assert_near(r.x, 3.0)
  runner.assert_near(r.y, 4.0)
  runner.assert_near(r.z, 0.0)

  local refl = vec3.reflect(v, n)
  runner.assert_near(refl.x, 3.0)
  runner.assert_near(refl.y, 4.0)
  runner.assert_near(refl.z, -5.0)
end)

runner.test("Skew-symmetric cross product matrix", function()
  local v = vec3.new(2.0, -3.0, 5.0)
  local w = vec3.new(1.0, 4.0, -2.0)
  local S = vec3.skew(v)

  -- S * w should equal v x w
  local Sw = vec3.new(
    S[1][1] * w.x + S[1][2] * w.y + S[1][3] * w.z,
    S[2][1] * w.x + S[2][2] * w.y + S[2][3] * w.z,
    S[3][1] * w.x + S[3][2] * w.y + S[3][3] * w.z
  )
  local v_cross_w = vec3.cross(v, w)
  runner.assert_near(Sw.x, v_cross_w.x)
  runner.assert_near(Sw.y, v_cross_w.y)
  runner.assert_near(Sw.z, v_cross_w.z)
end)
