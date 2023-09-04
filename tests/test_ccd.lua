local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local shapes = require("fulcrum.geometry.shapes")
local ccd = require("fulcrum.geometry.ccd")

runner.suite("Continuous Collision Detection (CCD)")

runner.test("Tunneling prevention and exact time of impact (TOI)", function()
  -- Sphere 1 static at origin, radius 0.5
  local s1 = shapes.sphere(0.5, vec3.new(0, 0, 0))
  local v1 = vec3.zero()

  -- Sphere 2 starting at (5, 0, 0), moving fast toward -X at 10 m/s, radius 0.5
  -- In discrete step dt=1.0, sphere 2 jumps to (-5, 0, 0), tunneling through s1
  local s2 = shapes.sphere(0.5, vec3.new(5, 0, 0))
  local v2 = vec3.new(-10, 0, 0)

  local hit, toi, normal, pt = ccd.solve_linear(s1, v1, s2, v2, 1.0, 1e-4)

  runner.assert_true(hit, "CCD detects high-speed collision")
  -- Analytical TOI: (5.0 - (0.5 + 0.5)) / 10.0 = 4.0 / 10.0 = 0.40 s
  runner.assert_near(toi, 0.40, 1e-3, "Time of impact matches analytical solution (0.4s)")
  runner.assert_true(pt ~= nil, "Contact point computed")
end)

runner.test("Missed trajectories with lateral offset", function()
  local s1 = shapes.sphere(0.5, vec3.new(0, 0, 0))
  local v1 = vec3.zero()

  -- Sphere 2 has lateral offset y=2.0 (distance > 0.5 + 0.5)
  local s2 = shapes.sphere(0.5, vec3.new(5, 2, 0))
  local v2 = vec3.new(-10, 0, 0)

  local hit, toi = ccd.solve_linear(s1, v1, s2, v2, 1.0, 1e-4)
  runner.assert_false(hit, "Separated trajectory does not report collision")
  runner.assert_equal(toi, 1.0, "Time advances to end of window")
end)

runner.test("General trajectory parameterized solver", function()
  -- Parametric trajectory
  local function traj_a(t)
    return shapes.sphere(0.5, vec3.new(0, 0, 0))
  end

  local function traj_b(t)
    -- Moving along parabola x(t) = 3 - 4*t
    return shapes.sphere(0.5, vec3.new(3.0 - 4.0 * t, 0, 0))
  end

  local hit, toi = ccd.solve_general(traj_a, traj_b, 4.0, 1.0, 1e-4)
  runner.assert_true(hit, "General trajectory solver detects impact")
  -- Distance from center is 3.0 - 4.0*t = 1.0 => 4.0*t = 2.0 => t = 0.5s
  runner.assert_near(toi, 0.5, 1e-3, "General solver TOI matches analytical 0.5s")
end)
