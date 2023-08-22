local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local shapes = require("fulcrum.geometry.shapes")
local gjk = require("fulcrum.geometry.gjk")

runner.suite("Gilbert-Johnson-Keerthi (GJK) Collision Detection")

runner.test("Sphere-sphere collision detection", function()
  local s1 = shapes.sphere(1.0, vec3.new(0.0, 0.0, 0.0))
  -- Touching/overlapping sphere at distance 1.5 < 2.0
  local s2 = shapes.sphere(1.0, vec3.new(1.5, 0.0, 0.0))
  -- Disjoint sphere at distance 3.0 > 2.0
  local s3 = shapes.sphere(1.0, vec3.new(3.0, 0.0, 0.0))

  local hit12 = gjk.intersect(s1, s2)
  runner.assert_true(hit12, "Spheres at distance 1.5 must intersect")

  local hit13 = gjk.intersect(s1, s3)
  runner.assert_false(hit13, "Spheres at distance 3.0 must not intersect")
end)

runner.test("Rotated box and capsule intersection", function()
  local R = mat3.from_axis_angle(vec3.unit_z(), 0.75)
  local b = shapes.box(vec3.new(1.0, 1.0, 1.0), { R = R, p = vec3.new(0, 0, 0) })

  -- Capsule passing through rotated box
  local cap_intersecting = shapes.capsule(0.2, vec3.new(0.5, 0.5, -2.0), vec3.new(0.5, 0.5, 2.0))
  local hit_cap = gjk.intersect(b, cap_intersecting)
  runner.assert_true(hit_cap, "Capsule passing through box must intersect")

  -- Capsule well clear of box
  local cap_clear = shapes.capsule(0.2, vec3.new(5.0, 5.0, -1.0), vec3.new(5.0, 5.0, 1.0))
  local hit_clear = gjk.intersect(b, cap_clear)
  runner.assert_false(hit_clear, "Clear capsule must not intersect")
end)
