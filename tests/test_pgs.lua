local runner = require("tests.test_framework")
local matrix = require("fulcrum.math.matrix")
local friction = require("fulcrum.contact.friction")
local pgs = require("fulcrum.contact.pgs")
local lemke = require("fulcrum.contact.lemke")

runner.suite("Projected Gauss-Seidel and Friction Cones")

runner.test("Coulomb friction cone projection", function()
  local mu = 0.4
  -- Case 1: inside cone (fn = 10, ft = 2 <= 4)
  local fn1, ft1, ft2 = friction.project_coulomb_cone(10.0, 1.2, 1.6, mu)
  runner.assert_near(fn1, 10.0, 1e-10)
  runner.assert_near(ft1, 1.2, 1e-10)
  runner.assert_near(ft2, 1.6, 1e-10)

  -- Case 2: outside cone (fn = 5, ft = 6 > 2)
  local fn2, ft2_1, ft2_2 = friction.project_coulomb_cone(5.0, 3.6, 4.8, mu)
  runner.assert_near(fn2, 5.0, 1e-10)
  local ft_norm = math.sqrt(ft2_1 * ft2_1 + ft2_2 * ft2_2)
  runner.assert_near(ft_norm, mu * fn2, 1e-10)
end)

runner.test("Projected Gauss-Seidel matches Lemke solution", function()
  local M = matrix.from_table({
    { 2.0, 1.0 },
    { 1.0, 2.0 }
  })
  local q = { -5.0, -6.0 }

  local _, z_lemke, _, _ = lemke.solve(M, q)
  local z_pgs, _ = pgs.solve_box(M, q, { 0.0, 0.0 }, { 1e10, 1e10 }, 200, 1e-8, 1.1)

  runner.assert_near(z_pgs[1], z_lemke[1], 1e-5)
  runner.assert_near(z_pgs[2], z_lemke[2], 1e-5)
end)
