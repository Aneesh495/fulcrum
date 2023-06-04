local runner = require("tests.test_framework")
local matrix = require("fulcrum.math.matrix")

runner.suite("Dense Matrix Algebra")

runner.test("Matrix creation and submatrix extraction", function()
  local M = matrix.zeros(4, 5)
  runner.assert_equal(M.rows, 4)
  runner.assert_equal(M.cols, 5)
  runner.assert_equal(M[2][3], 0.0)

  local A = matrix.from_table({
    { 1, 2, 3 },
    { 4, 5, 6 },
    { 7, 8, 9 }
  })
  local sub = matrix.submatrix(A, 1, 2, 2, 3)
  runner.assert_equal(sub.rows, 2)
  runner.assert_equal(sub.cols, 2)
  runner.assert_equal(sub[1][1], 2)
  runner.assert_equal(sub[1][2], 3)
  runner.assert_equal(sub[2][1], 5)
  runner.assert_equal(sub[2][2], 6)
end)

runner.test("Rectangular matrix multiplication", function()
  local A = matrix.from_table({
    { 1, 2, 3 },
    { 4, 5, 6 }
  }) -- 2x3
  local B = matrix.from_table({
    { 7, 8 },
    { 9, 1 },
    { 2, 3 }
  }) -- 3x2

  local C = A * B -- 2x2
  runner.assert_equal(C.rows, 2)
  runner.assert_equal(C.cols, 2)

  -- C[1][1] = 1*7 + 2*9 + 3*2 = 7 + 18 + 6 = 31
  -- C[1][2] = 1*8 + 2*1 + 3*3 = 8 + 2 + 9 = 19
  -- C[2][1] = 4*7 + 5*9 + 6*2 = 28 + 45 + 12 = 85
  -- C[2][2] = 4*8 + 5*1 + 6*3 = 32 + 5 + 18 = 55
  runner.assert_near(C[1][1], 31.0)
  runner.assert_near(C[1][2], 19.0)
  runner.assert_near(C[2][1], 85.0)
  runner.assert_near(C[2][2], 55.0)
end)

runner.test("LU linear solver and inversion", function()
  -- Hilbert-like 4x4 test matrix
  local A = matrix.from_table({
    { 10.0, -1.0,  2.0,  0.0 },
    { -1.0, 11.0, -1.0,  3.0 },
    {  2.0, -1.0, 10.0, -1.0 },
    {  0.0,  3.0, -1.0,  8.0 }
  })
  local b = { 6.0, 25.0, -11.0, 15.0 }

  local x = matrix.solve_lu(A, b)
  local Ax = matrix.mul_vec(A, x)
  runner.assert_vector_near(Ax, b, 1e-10)

  local Ainv = matrix.inv(A)
  local prod = A * Ainv
  local I = matrix.eye(4)

  for i = 1, 4 do
    for j = 1, 4 do
      runner.assert_near(prod[i][j], I[i][j], 1e-10)
    end
  end
end)

runner.test("Matrix determinant", function()
  local A = matrix.from_table({
    { 6, 1, 1 },
    { 4, -2, 5 },
    { 2, 8, 7 }
  })
  -- det(A) = 6*(-14 - 40) - 1*(28 - 10) + 1*(32 - (-4))
  --        = 6*(-54) - 18 + 36 = -324 - 18 + 36 = -306
  local d = matrix.det(A)
  runner.assert_near(d, -306.0, 1e-8)
end)
