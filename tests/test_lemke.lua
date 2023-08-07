local runner = require("tests.test_framework")
local matrix = require("fulcrum.math.matrix")
local lemke = require("fulcrum.contact.lemke")

runner.suite("Lemke Complementary Pivot LCP Solver")

runner.test("2D standard LCP test problem", function()
  -- M = [ 2, 1 ; 1, 2 ], q = [ -5, -6 ]
  local M = matrix.from_table({
    { 2.0, 1.0 },
    { 1.0, 2.0 }
  })
  local q = { -5.0, -6.0 }

  local status, z, w, iter = lemke.solve(M, q)
  runner.assert_equal(status, lemke.STATUS_SOLVED)

  -- Check w = M*z + q
  local Mz = matrix.mul_vec(M, z)
  for i = 1, 2 do
    runner.assert_near(w[i], Mz[i] + q[i], 1e-8)
    runner.assert_true(z[i] >= -1e-9, "z must be non-negative")
    runner.assert_true(w[i] >= -1e-9, "w must be non-negative")
    runner.assert_near(z[i] * w[i], 0.0, 1e-8) -- Complementarity
  end

  -- Exact analytical solution: z = [ 4/3, 7/3 ]
  runner.assert_near(z[1], 4.0 / 3.0, 1e-8)
  runner.assert_near(z[2], 7.0 / 3.0, 1e-8)
  runner.assert_near(w[1], 0.0, 1e-8)
  runner.assert_near(w[2], 0.0, 1e-8)
end)

runner.test("Partially active constraints with strictly positive slack", function()
  -- M = [ 2, 1 ; 1, 2 ], q = [ -2, 1 ]
  -- If z2 = 0, then 2*z1 - 2 = 0 -> z1 = 1.
  -- w1 = 0.
  -- w2 = 1*1 + 2*0 + 1 = 2 > 0 (inactive constraint).
  local M = matrix.from_table({
    { 2.0, 1.0 },
    { 1.0, 2.0 }
  })
  local q = { -2.0, 1.0 }

  local status, z, w, _ = lemke.solve(M, q)
  runner.assert_equal(status, lemke.STATUS_SOLVED)

  runner.assert_near(z[1], 1.0, 1e-8)
  runner.assert_near(z[2], 0.0, 1e-8)
  runner.assert_near(w[1], 0.0, 1e-8)
  runner.assert_near(w[2], 2.0, 1e-8)
end)

runner.test("Multi-contact 4D complementarity problem", function()
  local M = matrix.from_table({
    {  4, -1,  0,  0 },
    { -1,  4, -1,  0 },
    {  0, -1,  4, -1 },
    {  0,  0, -1,  4 }
  })
  local q = { -8.0, 2.0, -10.0, 4.0 }

  local status, z, w, _ = lemke.solve(M, q)
  runner.assert_equal(status, lemke.STATUS_SOLVED)

  local Mz = matrix.mul_vec(M, z)
  for i = 1, 4 do
    runner.assert_near(w[i], Mz[i] + q[i], 1e-8)
    runner.assert_true(z[i] >= -1e-9)
    runner.assert_true(w[i] >= -1e-9)
    runner.assert_near(z[i] * w[i], 0.0, 1e-8)
  end
end)
