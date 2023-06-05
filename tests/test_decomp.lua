local runner = require("tests.test_framework")
local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")

runner.suite("Matrix Decompositions")

runner.test("LU factorisation with partial pivoting", function()
  local A = matrix.from_table({
    { 2, 1, 1 },
    { 4, 3, 3 },
    { 8, 7, 9 }
  })
  local L, U, P = decomp.lu(A)

  -- Verify P * A = L * U
  local PA = P * A
  local LU = L * U

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(PA[i][j], LU[i][j], 1e-10)
    end
  end
end)

runner.test("Cholesky decomposition of SPD matrix", function()
  -- Symmetric positive-definite matrix
  local A = matrix.from_table({
    {  4,  12, -16 },
    { 12,  37, -43 },
    { -16, -43,  98 }
  })
  local L = decomp.cholesky(A)
  local Lt = matrix.transpose(L)
  local LLt = L * Lt

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(LLt[i][j], A[i][j], 1e-10)
    end
  end

  local b = { 8.0, 15.0, 23.0 }
  local x = decomp.cholesky_solve(L, b)
  local Ax = matrix.mul_vec(A, x)
  runner.assert_vector_near(Ax, b, 1e-10)
end)

runner.test("QR Householder decomposition", function()
  local A = matrix.from_table({
    { 12, -51,   4 },
    {  6, 167, -68 },
    { -4,  24, -41 }
  })
  local Q, R = decomp.qr(A)

  -- Check Q^T * Q = I
  local QtQ = matrix.transpose(Q) * Q
  local I = matrix.eye(3)
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(QtQ[i][j], I[i][j], 1e-10)
    end
  end

  -- Check Q * R = A
  local QR = Q * R
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(QR[i][j], A[i][j], 1e-10)
    end
  end
end)

runner.test("SVD and pseudo-inverse", function()
  local A = matrix.from_table({
    { 1, 2 },
    { 3, 4 },
    { 5, 6 }
  }) -- 3x2
  local U, S, V = decomp.svd(A)

  -- Reconstruct A_rec = U * diag(S) * V^T
  local m, n = A.rows, A.cols
  for i = 1, m do
    for j = 1, n do
      local sum = 0.0
      for k = 1, #S do
        sum = sum + U[i][k] * S[k] * V[j][k]
      end
      runner.assert_near(sum, A[i][j], 1e-7)
    end
  end

  -- Check pseudo-inverse condition: A * A^+ * A = A
  local Apinv = decomp.pinv(A)
  local rec = A * Apinv * A
  for i = 1, m do
    for j = 1, n do
      runner.assert_near(rec[i][j], A[i][j], 1e-7)
    end
  end
end)
