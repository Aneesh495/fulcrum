local runner = require("tests.test_framework")
local dual = require("fulcrum.math.dual")
local autodiff = require("fulcrum.math.autodiff")

runner.suite("Automatic Differentiation")

runner.test("Dual number arithmetic and elementary functions", function()
  -- f(x) = x^3 - 4*x^2 + 7*x - 5
  -- f'(x) = 3*x^2 - 8*x + 7
  local x_val = 3.0
  local x = dual.var(x_val, 1.0)
  local y = (x ^ 3.0) - (x ^ 2.0) * 4.0 + (x * 7.0) - 5.0

  local expected_val = 3^3 - 4*(3^2) + 7*3 - 5
  local expected_deriv = 3*(3^2) - 8*3 + 7

  runner.assert_near(y.r, expected_val, 1e-10)
  runner.assert_near(y.e, expected_deriv, 1e-10)

  -- Trigonometric derivative: d/dx [sin(2x)] = 2 * cos(2x)
  local z = dual.sin(x * 2.0)
  runner.assert_near(z.r, math.sin(6.0), 1e-10)
  runner.assert_near(z.e, 2.0 * math.cos(6.0), 1e-10)
end)

runner.test("Multivariable gradient of Rosenbrock function", function()
  -- f(x, y) = (1 - x)^2 + 100 * (y - x^2)^2
  -- df/dx = -2*(1 - x) - 400*x*(y - x^2)
  -- df/dy = 200*(y - x^2)
  local rosenbrock = function(v)
    local x = v[1]
    local y = v[2]
    local t1 = 1.0 - x
    local t2 = y - (x * x)
    return (t1 * t1) + (t2 * t2) * 100.0
  end

  local x0 = { 1.2, 2.5 }
  local grad = autodiff.grad(rosenbrock, x0)

  local x, y = x0[1], x0[2]
  local expected_df_dx = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x)
  local expected_df_dy = 200.0 * (y - x * x)

  runner.assert_near(grad[1], expected_df_dx, 1e-9)
  runner.assert_near(grad[2], expected_df_dy, 1e-9)
end)

runner.test("Vector Jacobian computation", function()
  -- f(x, y) = [ x^2 * y,  sin(x) * cos(y) ]
  -- J = [ 2*x*y,      x^2            ]
  --     [ cos(x)*cos(y), -sin(x)*sin(y) ]
  local vec_fn = function(v)
    local x, y = v[1], v[2]
    return {
      (x * x) * y,
      dual.sin(x) * dual.cos(y)
    }
  end

  local x0 = { 0.8, -1.5 }
  local J = autodiff.jacobian(vec_fn, x0)

  local x, y = x0[1], x0[2]
  runner.assert_near(J[1][1], 2.0 * x * y, 1e-10)
  runner.assert_near(J[1][2], x * x, 1e-10)
  runner.assert_near(J[2][1], math.cos(x) * math.cos(y), 1e-10)
  runner.assert_near(J[2][2], -math.sin(x) * math.sin(y), 1e-10)
end)

runner.test("Exact Hessian computation via hyper-dual numbers", function()
  -- f(x, y) = x^2 * y + sin(x)
  -- d^2f / dx^2 = 2*y - sin(x)
  -- d^2f / dy^2 = 0
  -- d^2f / dx dy = 2*x
  local hd = autodiff.hyper_dual()
  local f = function(v)
    local x, y = v[1], v[2]
    return (x * x) * y + hd.sin(x)
  end

  local x0 = { 1.3, -0.7 }
  local H = autodiff.hessian(f, x0)

  local x, y = x0[1], x0[2]
  runner.assert_near(H[1][1], 2.0 * y - math.sin(x), 1e-10)
  runner.assert_near(H[1][2], 2.0 * x, 1e-10)
  runner.assert_near(H[2][1], 2.0 * x, 1e-10)
  runner.assert_near(H[2][2], 0.0, 1e-10)
end)
