local runner = require("tests.test_framework")
local polynomial = require("fulcrum.math.polynomial")
local quadrature = require("fulcrum.math.quadrature")

runner.suite("Polynomials and Numerical Quadrature")

runner.test("Quadratic and cubic root finding", function()
  -- (x - 2)(x + 3) = x^2 + x - 6 = 0
  local q_roots = polynomial.roots_quadratic(1, 1, -6)
  table.sort(q_roots)
  runner.assert_near(q_roots[1], -3.0, 1e-10)
  runner.assert_near(q_roots[2], 2.0, 1e-10)

  -- (x - 1)(x - 2)(x - 3) = x^3 - 6*x^2 + 11*x - 6 = 0
  local c_roots = polynomial.roots_cubic(1, -6, 11, -6)
  table.sort(c_roots)
  runner.assert_equal(#c_roots, 3)
  runner.assert_near(c_roots[1], 1.0, 1e-9)
  runner.assert_near(c_roots[2], 2.0, 1e-9)
  runner.assert_near(c_roots[3], 3.0, 1e-9)
end)

runner.test("Quintic spline boundary condition preservation", function()
  local p0, v0, a0 = 0.0, 0.0, 0.0
  local p1, v1, a1 = 10.0, 0.0, 0.0
  local T = 2.5

  local spline = polynomial.quintic_spline(p0, v0, a0, p1, v1, a1, T)

  runner.assert_near(spline.pos(0.0), p0, 1e-10)
  runner.assert_near(spline.vel(0.0), v0, 1e-10)
  runner.assert_near(spline.acc(0.0), a0, 1e-10)

  runner.assert_near(spline.pos(T), p1, 1e-10)
  runner.assert_near(spline.vel(T), v1, 1e-10)
  runner.assert_near(spline.acc(T), a1, 1e-10)
end)

runner.test("Gauss-Legendre and Adaptive Simpson integration", function()
  -- Integrate f(x) = sin(x) on [0, pi] -> exact = 2.0
  local f = function(x) return math.sin(x) end
  local gl_res = quadrature.gauss_legendre(f, 0.0, math.pi, 6)
  runner.assert_near(gl_res, 2.0, 1e-9)

  local as_res = quadrature.adaptive_simpson(f, 0.0, math.pi, 1e-10)
  runner.assert_near(as_res, 2.0, 1e-9)
end)
