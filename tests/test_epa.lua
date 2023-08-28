local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local shapes = require("fulcrum.geometry.shapes")
local gjk = require("fulcrum.geometry.gjk")
local epa = require("fulcrum.geometry.epa")

runner.suite("Expanding Polytope Algorithm (EPA)")

runner.test("Sphere-sphere penetration depth and normal", function()
  local s1 = shapes.sphere(1.0, vec3.new(0, 0, 0))
  local s2 = shapes.sphere(1.0, vec3.new(1.6, 0, 0))

  local hit, simplex = gjk.intersect(s1, s2)
  runner.assert_true(hit, "Spheres should collide")

  local success, depth, normal = epa.penetration(s1, s2, simplex)
  runner.assert_true(success, "EPA should succeed")

  -- Expected depth: (1.0 + 1.0) - 1.6 = 0.4
  runner.assert_near(depth, 0.4, 1e-3, "Penetration depth matches analytical value")
  -- Normal should point along x-axis
  runner.assert_near(math.abs(normal.x), 1.0, 1e-2, "Normal x component is unitary")
  runner.assert_near(normal.y, 0.0, 1e-2, "Normal y component is zero")
  runner.assert_near(normal.z, 0.0, 1e-2, "Normal z component is zero")
end)

runner.test("Box-sphere contact penetration", function()
  -- Box centered at origin of half-extents (1, 1, 1)
  local b = shapes.box(vec3.new(1, 1, 1), vec3.new(0, 0, 0))
  -- Sphere centered at (1.3, 0, 0) of radius 0.5
  -- Penetration along x: (1.0 + 0.5) - 1.3 = 0.2
  local s = shapes.sphere(0.5, vec3.new(1.3, 0, 0))

  local hit, simplex = gjk.intersect(b, s)
  runner.assert_true(hit, "Box and sphere should collide")

  local success, depth, normal = epa.penetration(b, s, simplex)
  runner.assert_true(success, "EPA should converge on box-sphere penetration")
  runner.assert_near(depth, 0.2, 1e-3, "Box-sphere depth is 0.2")
end)
