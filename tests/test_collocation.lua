local runner = require("tests.test_framework")
local collocation = require("fulcrum.opt.collocation")

runner.suite("Direct Collocation (Hermite-Simpson)")

runner.test("Hermite-Simpson defect evaluates zero on exact cubic polynomial", function()
  -- Dynamics: dx/dt = 3 * t^2 => exact solution x(t) = t^3
  local h = 0.2
  local t_k = 1.0
  local t_next = t_k + h

  local x_k = { t_k^3 }
  local x_next = { t_next^3 }
  local u_k = { t_k }
  local u_next = { t_next }

  local function dynamics(_, u)
    local t = u[1]
    return { 3.0 * t * t }
  end

  local defect = collocation.hermite_simpson_defect(x_k, u_k, x_next, u_next, dynamics, h)
  -- Hermite-Simpson is 4th-order accurate and exact for polynomials up to degree 3
  runner.assert_near(defect[1], 0.0, 1e-6, "Defect vanishes identically for cubic trajectory")
end)

runner.test("Minimum energy state transfer via direct collocation", function()
  -- Double integrator: dx/dt = [v, u]
  local prob = {
    n_x = 2,
    n_u = 1,
    segments = 6,
    duration = 1.0,
    x_start = { 0.0, 0.0 },
    x_goal = { 1.0, 0.0 },
    dynamics = function(x, u)
      local v = x[2]
      local a = u[1]
      return { v, a }
    end,
    cost = function(_, u)
      local a = u[1]
      return 0.5 * a * a
    end
  }

  local res = collocation.solve(prob, {
    max_iter = 10,
    tol_defect = 1e-4
  })

  -- Verify boundary conditions
  local x_end = res.x[#res.x]
  runner.assert_near(x_end[1], 1.0, 0.01, "Final position reaches target 1.0")
  runner.assert_near(x_end[2], 0.0, 0.01, "Final velocity near target 0.0")

  -- Verify defect constraints across segments
  for k = 1, prob.segments do
    local def = collocation.hermite_simpson_defect(
      res.x[k], res.u[k], res.x[k + 1], res.u[k + 1], prob.dynamics, res.dt
    )
    runner.assert_near(def[1], 0.0, 0.01, "Position defect within tolerance")
    runner.assert_near(def[2], 0.0, 0.01, "Velocity defect within tolerance")
  end
end)
