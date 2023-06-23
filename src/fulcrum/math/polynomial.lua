local polynomial = {}

function polynomial.eval(coeffs, x)
  -- Horner scheme: coeffs = { a_n, a_{n-1}, ..., a_0 } for a_n * x^n + ... + a_0
  local val = coeffs[1] or 0.0
  for i = 2, #coeffs do
    val = val * x + coeffs[i]
  end
  return val
end

function polynomial.derivative(coeffs)
  local n = #coeffs - 1
  if n <= 0 then return { 0.0 } end
  local d = {}
  for i = 1, n do
    d[i] = coeffs[i] * (n - i + 1)
  end
  return d
end

function polynomial.roots_quadratic(a, b, c)
  -- Solves a*x^2 + b*x + c = 0
  if math.abs(a) < 1e-15 then
    if math.abs(b) < 1e-15 then return {} end
    return { -c / b }
  end

  local disc = b * b - 4.0 * a * c
  if disc > 1e-12 then
    local sqrt_disc = math.sqrt(disc)
    local q = -0.5 * (b + (b >= 0 and sqrt_disc or -sqrt_disc))
    return { q / a, c / q }
  elseif math.abs(disc) <= 1e-12 then
    return { -0.5 * b / a }
  else
    -- Complex roots (return empty list of real roots)
    return {}
  end
end

function polynomial.roots_cubic(a, b, c, d)
  -- Cardano formula for real roots of a*x^3 + b*x^2 + c*x + d = 0
  if math.abs(a) < 1e-15 then
    return polynomial.roots_quadratic(b, c, d)
  end

  -- Depress cubic: t = x + b / (3*a) -> t^3 + p*t + q = 0
  local inv_a = 1.0 / a
  local b_a = b * inv_a
  local c_a = c * inv_a
  local d_a = d * inv_a

  local p = c_a - (b_a * b_a) / 3.0
  local q = (2.0 * b_a * b_a * b_a) / 27.0 - (b_a * c_a) / 3.0 + d_a
  local shift = -b_a / 3.0

  local disc = (q * q) * 0.25 + (p * p * p) / 27.0

  if disc > 1e-12 then
    -- One real root
    local sqrt_disc = math.sqrt(disc)
    local u = -0.5 * q + sqrt_disc
    local v = -0.5 * q - sqrt_disc
    local cbrt_u = (u >= 0) and (u ^ (1.0 / 3.0)) or (-((-u) ^ (1.0 / 3.0)))
    local cbrt_v = (v >= 0) and (v ^ (1.0 / 3.0)) or (-((-v) ^ (1.0 / 3.0)))
    return { cbrt_u + cbrt_v + shift }
  elseif math.abs(disc) <= 1e-12 then
    -- Multiple real roots
    local u = -0.5 * q
    local cbrt_u = (u >= 0) and (u ^ (1.0 / 3.0)) or (-((-u) ^ (1.0 / 3.0)))
    return { 2.0 * cbrt_u + shift, -cbrt_u + shift }
  else
    -- Three distinct real roots (casus irreducibilis via trigonometry)
    local r = math.sqrt(-(p * p * p) / 27.0)
    local phi = math.acos(math.max(-1.0, math.min(1.0, -0.5 * q / r)))
    local m = 2.0 * (r ^ (1.0 / 3.0))
    local r1 = m * math.cos(phi / 3.0) + shift
    local r2 = m * math.cos((phi + 2.0 * math.pi) / 3.0) + shift
    local r3 = m * math.cos((phi + 4.0 * math.pi) / 3.0) + shift
    return { r1, r2, r3 }
  end
end

-- Quintic Hermite spline for minimum-jerk trajectory generation
-- Specified by (p0, v0, a0) at t=0 and (p1, v1, a1) at t=T
function polynomial.quintic_spline(p0, v0, a0, p1, v1, a1, T)
  local T2 = T * T
  local T3 = T2 * T
  local T4 = T3 * T
  local T5 = T4 * T

  local h = p1 - p0
  local a0_term = a0
  local a1_term = a1

  local c0 = p0
  local c1 = v0
  local c2 = 0.5 * a0

  local c3 = (20.0 * h - (8.0 * v1 + 12.0 * v0) * T - (3.0 * a0_term - a1_term) * T2) / (2.0 * T3)
  local c4 = (-30.0 * h + (14.0 * v1 + 16.0 * v0) * T + (3.0 * a0_term - 2.0 * a1_term) * T2) / (2.0 * T4)
  local c5 = (12.0 * h - 6.0 * (v1 + v0) * T + (a1_term - a0_term) * T2) / (2.0 * T5)

  local spline = {}
  function spline.pos(t)
    local s = math.max(0.0, math.min(T, t))
    local s2 = s * s
    local s3 = s2 * s
    local s4 = s3 * s
    local s5 = s4 * s
    return c0 + c1 * s + c2 * s2 + c3 * s3 + c4 * s4 + c5 * s5
  end
  function spline.vel(t)
    local s = math.max(0.0, math.min(T, t))
    local s2 = s * s
    local s3 = s2 * s
    local s4 = s3 * s
    return c1 + 2.0 * c2 * s + 3.0 * c3 * s2 + 4.0 * c4 * s3 + 5.0 * c5 * s4
  end
  function spline.acc(t)
    local s = math.max(0.0, math.min(T, t))
    local s2 = s * s
    local s3 = s2 * s
    return 2.0 * c2 + 6.0 * c3 * s + 12.0 * c4 * s2 + 20.0 * c5 * s3
  end

  return spline
end

return polynomial
