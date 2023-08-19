local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local shapes = require("fulcrum.geometry.shapes")

runner.suite("Convex Shapes and Support Functions")

runner.test("Sphere support mapping", function()
  local s = shapes.sphere(2.5, vec3.new(1.0, -1.0, 3.0))
  local d = vec3.unit_x()
  local p = s.support(d)

  runner.assert_near(p.x, 3.5, 1e-10)
  runner.assert_near(p.y, -1.0, 1e-10)
  runner.assert_near(p.z, 3.0, 1e-10)
end)

runner.test("Transformed box support mapping", function()
  local h = vec3.new(1.0, 2.0, 3.0)
  local R = mat3.from_axis_angle(vec3.unit_z(), math.pi * 0.5) -- 90 deg around Z
  local p_trans = vec3.new(5.0, 10.0, 0.0)

  local b = shapes.box(h, { R = R, p = p_trans })
  -- In world frame, 90 deg around Z swaps X and Y:
  -- d = [1, 0, 0] in world -> d_local = R^T * d = [0, -1, 0]
  -- local box support for d_local is [hx, -hy, hz] = [1, -2, 3]
  -- World support = R * [1, -2, 3] + p_trans = [2, 1, 3] + [5, 10, 0] = [7, 11, 3]
  local p = b.support(vec3.unit_x())
  runner.assert_near(p.x, 7.0, 1e-10)
  runner.assert_near(p.y, 11.0, 1e-10)
  runner.assert_near(p.z, 3.0, 1e-10)
end)

runner.test("Ellipsoid analytical support function", function()
  local e = shapes.ellipsoid(vec3.new(2.0, 1.0, 1.0))
  local d = vec3.unit_x()
  local p = e.support(d)

  runner.assert_near(p.x, 2.0, 1e-10)
  runner.assert_near(p.y, 0.0, 1e-10)
  runner.assert_near(p.z, 0.0, 1e-10)
end)
