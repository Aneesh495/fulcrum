local matrix = {}
local matrix_mt = { __index = matrix }

function matrix.new(rows, cols, init)
  local m = {
    rows = rows,
    cols = cols
  }
  for i = 1, rows do
    m[i] = {}
    for j = 1, cols do
      if type(init) == "function" then
        m[i][j] = init(i, j)
      elseif type(init) == "number" then
        m[i][j] = init
      else
        m[i][j] = 0.0
      end
    end
  end
  return setmetatable(m, matrix_mt)
end

function matrix.zeros(rows, cols)
  return matrix.new(rows, cols, 0.0)
end

function matrix.ones(rows, cols)
  return matrix.new(rows, cols, 1.0)
end

function matrix.eye(n)
  local m = matrix.zeros(n, n)
  for i = 1, n do
    m[i][i] = 1.0
  end
  return m
end

function matrix.diag(v)
  local n = #v
  local m = matrix.zeros(n, n)
  for i = 1, n do
    m[i][i] = v[i]
  end
  return m
end

function matrix.from_table(tbl)
  local rows = #tbl
  local cols = #tbl[1]
  local m = matrix.new(rows, cols)
  for i = 1, rows do
    for j = 1, cols do
      m[i][j] = tbl[i][j] or 0.0
    end
  end
  return m
end

function matrix.clone(m)
  local res = matrix.new(m.rows, m.cols)
  for i = 1, m.rows do
    for j = 1, m.cols do
      res[i][j] = m[i][j]
    end
  end
  return res
end

function matrix.transpose(m)
  local res = matrix.new(m.cols, m.rows)
  for i = 1, m.rows do
    for j = 1, m.cols do
      res[j][i] = m[i][j]
    end
  end
  return res
end

function matrix.trace(m)
  local n = math.min(m.rows, m.cols)
  local tr = 0.0
  for i = 1, n do
    tr = tr + m[i][i]
  end
  return tr
end

function matrix.submatrix(m, r1, r2, c1, c2)
  local rows = r2 - r1 + 1
  local cols = c2 - c1 + 1
  local res = matrix.new(rows, cols)
  for i = 1, rows do
    for j = 1, cols do
      res[i][j] = m[r1 + i - 1][c1 + j - 1]
    end
  end
  return res
end

function matrix.set_submatrix(dst, r, c, src)
  for i = 1, src.rows do
    for j = 1, src.cols do
      dst[r + i - 1][c + j - 1] = src[i][j]
    end
  end
  return dst
end

function matrix.mul(A, B)
  if A.cols ~= B.rows then
    error(string.format("Matrix dimension mismatch in multiplication: (%dx%d) * (%dx%d)",
      A.rows, A.cols, B.rows, B.cols))
  end
  local res = matrix.zeros(A.rows, B.cols)
  for i = 1, A.rows do
    local ai = A[i]
    local resi = res[i]
    for k = 1, A.cols do
      local aik = ai[k]
      if aik ~= 0.0 then
        local bk = B[k]
        for j = 1, B.cols do
          resi[j] = resi[j] + aik * bk[j]
        end
      end
    end
  end
  return res
end

function matrix.mul_vec(A, x)
  if A.cols ~= #x then
    error(string.format("Matrix-vector dimension mismatch: (%dx%d) * (%d)",
      A.rows, A.cols, #x))
  end
  local res = {}
  for i = 1, A.rows do
    local sum = 0.0
    local ai = A[i]
    for j = 1, A.cols do
      sum = sum + ai[j] * x[j]
    end
    res[i] = sum
  end
  return res
end

function matrix.solve_lu(A, b)
  if A.rows ~= A.cols then
    error("Matrix must be square to solve via LU decomposition")
  end
  local n = A.rows
  if #b ~= n then
    error("Right-hand side dimension must match matrix rows")
  end

  -- Clone matrix for in-place LU factorisation with partial pivoting
  local LU = matrix.clone(A)
  local piv = {}
  for i = 1, n do piv[i] = i end

  for k = 1, n do
    -- Find pivot
    local max_val = 0.0
    local p = k
    for i = k, n do
      local val = math.abs(LU[i][k])
      if val > max_val then
        max_val = val
        p = i
      end
    end

    if max_val < 1e-15 then
      error("Matrix is singular to machine precision in LU solve")
    end

    -- Swap rows if needed
    if p ~= k then
      piv[k], piv[p] = piv[p], piv[k]
      LU[k], LU[p] = LU[p], LU[k]
    end

    local piv_inv = 1.0 / LU[k][k]
    for i = k + 1, n do
      local factor = LU[i][k] * piv_inv
      LU[i][k] = factor
      for j = k + 1, n do
        LU[i][j] = LU[i][j] - factor * LU[k][j]
      end
    end
  end

  -- Forward substitution with permutation vector
  local y = {}
  for i = 1, n do
    local sum = b[piv[i]]
    for j = 1, i - 1 do
      sum = sum - LU[i][j] * y[j]
    end
    y[i] = sum
  end

  -- Backward substitution
  local x = {}
  for i = n, 1, -1 do
    local sum = y[i]
    for j = i + 1, n do
      sum = sum - LU[i][j] * x[j]
    end
    x[i] = sum / LU[i][i]
  end

  return x
end

function matrix.inv(A)
  if A.rows ~= A.cols then
    error("Only square matrices have an inverse")
  end
  local n = A.rows
  local res = matrix.zeros(n, n)

  for j = 1, n do
    local e = {}
    for i = 1, n do e[i] = (i == j) and 1.0 or 0.0 end
    local col = matrix.solve_lu(A, e)
    for i = 1, n do
      res[i][j] = col[i]
    end
  end

  return res
end

function matrix.det(A)
  if A.rows ~= A.cols then
    error("Determinant is defined only for square matrices")
  end
  local n = A.rows
  local LU = matrix.clone(A)
  local sign = 1.0

  for k = 1, n do
    local max_val = 0.0
    local p = k
    for i = k, n do
      local val = math.abs(LU[i][k])
      if val > max_val then
        max_val = val
        p = i
      end
    end

    if max_val < 1e-15 then
      return 0.0
    end

    if p ~= k then
      LU[k], LU[p] = LU[p], LU[k]
      sign = -sign
    end

    local piv_inv = 1.0 / LU[k][k]
    for i = k + 1, n do
      local factor = LU[i][k] * piv_inv
      LU[i][k] = factor
      for j = k + 1, n do
        LU[i][j] = LU[i][j] - factor * LU[k][j]
      end
    end
  end

  local d = sign
  for i = 1, n do
    d = d * LU[i][i]
  end
  return d
end

function matrix.norm_frobenius(A)
  local sum = 0.0
  for i = 1, A.rows do
    for j = 1, A.cols do
      local v = A[i][j]
      sum = sum + v * v
    end
  end
  return math.sqrt(sum)
end

-- Metamethods
function matrix_mt.__add(A, B)
  if A.rows ~= B.rows or A.cols ~= B.cols then
    error("Dimension mismatch in matrix addition")
  end
  local res = matrix.new(A.rows, A.cols)
  for i = 1, A.rows do
    for j = 1, A.cols do
      res[i][j] = A[i][j] + B[i][j]
    end
  end
  return res
end

function matrix_mt.__sub(A, B)
  if A.rows ~= B.rows or A.cols ~= B.cols then
    error("Dimension mismatch in matrix subtraction")
  end
  local res = matrix.new(A.rows, A.cols)
  for i = 1, A.rows do
    for j = 1, A.cols do
      res[i][j] = A[i][j] - B[i][j]
    end
  end
  return res
end

function matrix_mt.__mul(A, B)
  if type(A) == "number" then
    local res = matrix.new(B.rows, B.cols)
    for i = 1, B.rows do
      for j = 1, B.cols do
        res[i][j] = B[i][j] * A
      end
    end
    return res
  elseif type(B) == "number" then
    local res = matrix.new(A.rows, A.cols)
    for i = 1, A.rows do
      for j = 1, A.cols do
        res[i][j] = A[i][j] * B
      end
    end
    return res
  elseif type(B) == "table" and not B.rows then
    -- Matrix-vector multiplication
    return matrix.mul_vec(A, B)
  else
    return matrix.mul(A, B)
  end
end

function matrix_mt.__unm(A)
  local res = matrix.new(A.rows, A.cols)
  for i = 1, A.rows do
    for j = 1, A.cols do
      res[i][j] = -A[i][j]
    end
  end
  return res
end

function matrix_mt.__tostring(m)
  local lines = {}
  table.insert(lines, string.format("matrix(%dx%d):", m.rows, m.cols))
  for i = 1, m.rows do
    local row_str = {}
    for j = 1, m.cols do
      table.insert(row_str, string.format("%8.4f", m[i][j]))
    end
    table.insert(lines, "  [" .. table.concat(row_str, ", ") .. "]")
  end
  return table.concat(lines, "\n")
end

return matrix
