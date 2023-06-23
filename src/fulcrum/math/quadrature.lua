local quadrature = {}

-- Gauss-Legendre quadrature points and weights on [-1, 1]
local gl_data = {
  [2] = {
    x = { -0.5773502691896257, 0.5773502691896257 },
    w = { 1.0, 1.0 }
  },
  [4] = {
    x = {
      -0.8611363115940526, -0.3399810435848563,
       0.3399810435848563,  0.8611363115940526
    },
    w = {
      0.3478548451374538, 0.6521451548625461,
      0.6521451548625461, 0.3478548451374538
    }
  },
  [6] = {
    x = {
      -0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
       0.2386191860831969,  0.6612093864662645,  0.9324695142031521
    },
    w = {
      0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
      0.4679139345726910, 0.3607615730481386, 0.1713244923791704
    }
  }
}

function quadrature.gauss_legendre(f, a, b, order)
  order = order or 4
  local data = gl_data[order] or gl_data[4]
  local half_len = 0.5 * (b - a)
  local mid = 0.5 * (b + a)

  local sum = 0.0
  for i = 1, #data.x do
    local xi = mid + half_len * data.x[i]
    sum = sum + data.w[i] * f(xi)
  end

  return sum * half_len
end

function quadrature.adaptive_simpson(f, a, b, tol, max_depth)
  tol = tol or 1e-8
  max_depth = max_depth or 20

  local function simpson(left, right, f_left, f_mid, f_right)
    return (right - left) * (f_left + 4.0 * f_mid + f_right) / 6.0
  end

  local function step(left, right, eps, s, fa, fb, fc, depth)
    local mid = (left + right) * 0.5
    local d = (left + mid) * 0.5
    local e = (mid + right) * 0.5
    local fd = f(d)
    local fe = f(e)

    local sl = simpson(left, mid, fa, fd, fc)
    local sr = simpson(mid, right, fc, fe, fb)
    local s2 = sl + sr
    local err = s2 - s

    if depth >= max_depth or math.abs(err) <= 15.0 * eps then
      return s2 + err / 15.0
    end

    return step(left, mid, eps * 0.5, sl, fa, fc, fd, depth + 1) +
           step(mid, right, eps * 0.5, sr, fc, fb, fe, depth + 1)
  end

  local fa = f(a)
  local fb = f(b)
  local c = (a + b) * 0.5
  local fc = f(c)
  local s0 = simpson(a, b, fa, fc, fb)

  return step(a, b, tol, s0, fa, fb, fc, 1)
end

return quadrature
