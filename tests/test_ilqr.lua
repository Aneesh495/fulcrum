local runner = require("tests.test_framework")
local ilqr = require("fulcrum.opt.ilqr")

runner.suite("Iterative Linear Quadratic Regulator (iLQR / DDP)")

runner.test("Double integrator optimal regulation", function()
  -- State: x = [position, velocity]
  -- Control: u = [force/acceleration]
  local dt = 0.1
  local N = 25

  local prob = {
    n_x = 2,
    n_u = 1,
    N = N,
    dt = dt,
    dynamics = function(x, u, _)
      local p, v = x[1], x[2]
      local a = u[1]
      return {
        p + v * dt + 0.5 * a * dt * dt,
        v + a * dt
      }
    end,
    cost = function(x, u, _)
      -- Stage cost: penalize state deviation and control effort
      local p, v = x[1], x[2]
      local a = u[1]
      return 0.5 * (10.0 * p * p + 1.0 * v * v + 0.1 * a * a)
    end,
    terminal_cost = function(x)
      local p, v = x[1], x[2]
      return 0.5 * (50.0 * p * p + 10.0 * v * v)
    end
  }

  local x0 = { 2.0, 0.0 }
  local res = ilqr.solve(prob, x0, nil, {
    max_iter = 30,
    tol = 1e-4
  })

  runner.assert_true(#res.cost_history >= 2, "Performed at least one optimization iteration")
  local initial_cost = res.cost_history[1]
  local final_cost = res.cost
  runner.assert_true(final_cost < initial_cost * 0.25, "iLQR reduced total cost by over 75%")

  -- Final state should be close to origin
  local x_final = res.x[N + 1]
  runner.assert_near(x_final[1], 0.0, 0.15, "Final position driven to origin")
  runner.assert_near(x_final[2], 0.0, 0.15, "Final velocity driven to zero")
end)

runner.test("Nonlinear pendulum stabilization", function()
  -- State: x = [theta, omega]
  local dt = 0.05
  local N = 30
  local g = 9.81
  local l = 1.0

  local prob = {
    n_x = 2,
    n_u = 1,
    N = N,
    dt = dt,
    dynamics = function(x, u, _)
      local th, w = x[1], x[2]
      local tau = u[1]
      local w_dot = -(g / l) * math.sin(th) + tau
      return {
        th + w * dt,
        w + w_dot * dt
      }
    end,
    cost = function(x, u, _)
      local th, w = x[1], x[2]
      local tau = u[1]
      return 0.5 * (20.0 * th * th + 2.0 * w * w + 0.05 * tau * tau)
    end,
    terminal_cost = function(x)
      local th, w = x[1], x[2]
      return 0.5 * (100.0 * th * th + 20.0 * w * w)
    end
  }

  -- Perturbed starting angle (0.4 rad)
  local x0 = { 0.4, 0.0 }
  local res = ilqr.solve(prob, x0, nil, {
    max_iter = 25,
    tol = 1e-3
  })

  local initial_cost = res.cost_history[1]
  local final_cost = res.cost
  runner.assert_true(final_cost < initial_cost * 0.2, "Pendulum cost reduced significantly")
  local x_final = res.x[N + 1]
  runner.assert_near(x_final[1], 0.0, 0.1, "Pendulum angle stabilized near zero")
end)
