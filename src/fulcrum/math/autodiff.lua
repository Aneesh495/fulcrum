local dual = require("fulcrum.math.dual")
local matrix = require("fulcrum.math.matrix")

local autodiff = {}

function autodiff.grad(f, x)
  local n = #x
  local grad = {}

  for i = 1, n do
    local dual_x = {}
    for j = 1, n do
      dual_x[j] = dual.new(x[j], (i == j) and 1.0 or 0.0)
    end
    local res = f(dual_x)
    grad[i] = res.e
  end

  return grad
end

function autodiff.jacobian(f, x)
  local n = #x
  local m = nil
  local J = nil

  for j = 1, n do
    local dual_x = {}
    for k = 1, n do
      dual_x[k] = dual.new(x[k], (j == k) and 1.0 or 0.0)
    end
    local res = f(dual_x)
    if not J then
      m = #res
      J = matrix.zeros(m, n)
    end
    for i = 1, m do
      J[i][j] = res[i].e
    end
  end

  return J
end

function autodiff.directional(f, x, v)
  local n = #x
  local dual_x = {}
  for i = 1, n do
    dual_x[i] = dual.new(x[i], v[i])
  end
  local res = f(dual_x)
  return res.e
end

-- Hyper-dual number implementation for exact Hessian computation
local hyper_dual = {}
local hyper_dual_mt = { __index = hyper_dual }

function hyper_dual.new(r, e1, e2, e12)
  local h = {
    r = r or 0.0,
    e1 = e1 or 0.0,
    e2 = e2 or 0.0,
    e12 = e12 or 0.0
  }
  return setmetatable(h, hyper_dual_mt)
end

function hyper_dual_mt.__add(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae1 = type(a) == "table" and a.e1 or 0.0
  local ae2 = type(a) == "table" and a.e2 or 0.0
  local ae12 = type(a) == "table" and a.e12 or 0.0

  local br = type(b) == "table" and b.r or b
  local be1 = type(b) == "table" and b.e1 or 0.0
  local be2 = type(b) == "table" and b.e2 or 0.0
  local be12 = type(b) == "table" and b.e12 or 0.0

  return hyper_dual.new(ar + br, ae1 + be1, ae2 + be2, ae12 + be12)
end

function hyper_dual_mt.__sub(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae1 = type(a) == "table" and a.e1 or 0.0
  local ae2 = type(a) == "table" and a.e2 or 0.0
  local ae12 = type(a) == "table" and a.e12 or 0.0

  local br = type(b) == "table" and b.r or b
  local be1 = type(b) == "table" and b.e1 or 0.0
  local be2 = type(b) == "table" and b.e2 or 0.0
  local be12 = type(b) == "table" and b.e12 or 0.0

  return hyper_dual.new(ar - br, ae1 - be1, ae2 - be2, ae12 - be12)
end

function hyper_dual_mt.__mul(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae1 = type(a) == "table" and a.e1 or 0.0
  local ae2 = type(a) == "table" and a.e2 or 0.0
  local ae12 = type(a) == "table" and a.e12 or 0.0

  local br = type(b) == "table" and b.r or b
  local be1 = type(b) == "table" and b.e1 or 0.0
  local be2 = type(b) == "table" and b.e2 or 0.0
  local be12 = type(b) == "table" and b.e12 or 0.0

  return hyper_dual.new(
    ar * br,
    ar * be1 + ae1 * br,
    ar * be2 + ae2 * br,
    ar * be12 + ae1 * be2 + ae2 * be1 + ae12 * br
  )
end

function hyper_dual.sin(a)
  if type(a) == "number" then return math.sin(a) end
  local sa = math.sin(a.r)
  local ca = math.cos(a.r)
  return hyper_dual.new(
    sa,
    a.e1 * ca,
    a.e2 * ca,
    a.e12 * ca - a.e1 * a.e2 * sa
  )
end

function hyper_dual.cos(a)
  if type(a) == "number" then return math.cos(a) end
  local sa = math.sin(a.r)
  local ca = math.cos(a.r)
  return hyper_dual.new(
    ca,
    -a.e1 * sa,
    -a.e2 * sa,
    -a.e12 * sa - a.e1 * a.e2 * ca
  )
end

function hyper_dual.exp(a)
  if type(a) == "number" then return math.exp(a) end
  local ea = math.exp(a.r)
  return hyper_dual.new(
    ea,
    a.e1 * ea,
    a.e2 * ea,
    a.e12 * ea + a.e1 * a.e2 * ea
  )
end

function autodiff.hyper_dual()
  return hyper_dual
end

function autodiff.hessian(f, x)
  local n = #x
  local H = matrix.zeros(n, n)

  for i = 1, n do
    for j = i, n do
      local hx = {}
      for k = 1, n do
        local e1 = (k == i) and 1.0 or 0.0
        local e2 = (k == j) and 1.0 or 0.0
        hx[k] = hyper_dual.new(x[k], e1, e2, 0.0)
      end
      local res = f(hx)
      H[i][j] = res.e12
      H[j][i] = res.e12
    end
  end

  return H
end

return autodiff
