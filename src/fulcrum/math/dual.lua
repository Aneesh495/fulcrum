local dual = {}
local dual_mt = { __index = dual }

function dual.new(real, eps)
  local d = {
    r = real or 0.0,
    e = eps or 0.0
  }
  return setmetatable(d, dual_mt)
end

function dual.const(c)
  return dual.new(c, 0.0)
end

function dual.var(val, deriv)
  return dual.new(val, deriv or 1.0)
end

-- Math functions
function dual.sin(a)
  if type(a) == "number" then return math.sin(a) end
  return dual.new(math.sin(a.r), a.e * math.cos(a.r))
end

function dual.cos(a)
  if type(a) == "number" then return math.cos(a) end
  return dual.new(math.cos(a.r), -a.e * math.sin(a.r))
end

function dual.tan(a)
  if type(a) == "number" then return math.tan(a) end
  local cos_r = math.cos(a.r)
  return dual.new(math.tan(a.r), a.e / (cos_r * cos_r))
end

function dual.exp(a)
  if type(a) == "number" then return math.exp(a) end
  local exp_r = math.exp(a.r)
  return dual.new(exp_r, a.e * exp_r)
end

function dual.log(a)
  if type(a) == "number" then return math.log(a) end
  return dual.new(math.log(a.r), a.e / a.r)
end

function dual.sqrt(a)
  if type(a) == "number" then return math.sqrt(a) end
  local sqrt_r = math.sqrt(a.r)
  return dual.new(sqrt_r, 0.5 * a.e / sqrt_r)
end

function dual.abs(a)
  if type(a) == "number" then return math.abs(a) end
  local s = (a.r >= 0) and 1.0 or -1.0
  return dual.new(math.abs(a.r), a.e * s)
end

function dual.asin(a)
  if type(a) == "number" then return math.asin(a) end
  return dual.new(math.asin(a.r), a.e / math.sqrt(1.0 - a.r * a.r))
end

function dual.acos(a)
  if type(a) == "number" then return math.acos(a) end
  return dual.new(math.acos(a.r), -a.e / math.sqrt(1.0 - a.r * a.r))
end

function dual.atan(a)
  if type(a) == "number" then return math.atan(a) end
  return dual.new(math.atan(a.r), a.e / (1.0 + a.r * a.r))
end

function dual.atan2(y, x)
  local y_r = type(y) == "table" and y.r or y
  local y_e = type(y) == "table" and y.e or 0.0
  local x_r = type(x) == "table" and x.r or x
  local x_e = type(x) == "table" and x.e or 0.0

  local denom = x_r * x_r + y_r * y_r
  local deriv = (x_r * y_e - y_r * x_e) / denom
  return dual.new(math.atan(y_r, x_r), deriv)
end

-- Metamethods
function dual_mt.__add(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae = type(a) == "table" and a.e or 0.0
  local br = type(b) == "table" and b.r or b
  local be = type(b) == "table" and b.e or 0.0
  return dual.new(ar + br, ae + be)
end

function dual_mt.__sub(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae = type(a) == "table" and a.e or 0.0
  local br = type(b) == "table" and b.r or b
  local be = type(b) == "table" and b.e or 0.0
  return dual.new(ar - br, ae - be)
end

function dual_mt.__mul(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae = type(a) == "table" and a.e or 0.0
  local br = type(b) == "table" and b.r or b
  local be = type(b) == "table" and b.e or 0.0
  -- Product rule: (ar + ae eps)(br + be eps) = ar*br + (ar*be + ae*br) eps
  return dual.new(ar * br, ar * be + ae * br)
end

function dual_mt.__div(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae = type(a) == "table" and a.e or 0.0
  local br = type(b) == "table" and b.r or b
  local be = type(b) == "table" and b.e or 0.0
  -- Quotient rule: (ar + ae eps) / (br + be eps) = ar/br + ((ae*br - ar*be) / br^2) eps
  local inv_br = 1.0 / br
  return dual.new(ar * inv_br, (ae * br - ar * be) * (inv_br * inv_br))
end

function dual_mt.__pow(a, b)
  local ar = type(a) == "table" and a.r or a
  local ae = type(a) == "table" and a.e or 0.0
  local br = type(b) == "table" and b.r or b
  local be = type(b) == "table" and b.e or 0.0

  local pow_val = ar ^ br
  local deriv = 0.0
  if ae ~= 0.0 then
    deriv = deriv + br * (ar ^ (br - 1.0)) * ae
  end
  if be ~= 0.0 and ar > 0.0 then
    deriv = deriv + pow_val * math.log(ar) * be
  end
  return dual.new(pow_val, deriv)
end

function dual_mt.__unm(a)
  return dual.new(-a.r, -a.e)
end

function dual_mt.__tostring(a)
  return string.format("dual(%.6f, %.6f)", a.r, a.e)
end

return dual
