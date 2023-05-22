local vec3 = {}
local vec3_mt = { __index = vec3 }

function vec3.new(x, y, z)
  local v = {
    x = x or 0.0,
    y = y or 0.0,
    z = z or 0.0
  }
  return setmetatable(v, vec3_mt)
end

function vec3.zero()
  return vec3.new(0.0, 0.0, 0.0)
end

function vec3.unit_x()
  return vec3.new(1.0, 0.0, 0.0)
end

function vec3.unit_y()
  return vec3.new(0.0, 1.0, 0.0)
end

function vec3.unit_z()
  return vec3.new(0.0, 0.0, 1.0)
end

function vec3.clone(v)
  return vec3.new(v.x, v.y, v.z)
end

function vec3.set(v, x, y, z)
  v.x = x
  v.y = y
  v.z = z
  return v
end

function vec3.copy(dst, src)
  dst.x = src.x
  dst.y = src.y
  dst.z = src.z
  return dst
end

function vec3.unpack(v)
  return v.x, v.y, v.z
end

function vec3.dot(a, b)
  return a.x * b.x + a.y * b.y + a.z * b.z
end

function vec3.cross(a, b)
  return vec3.new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  )
end

function vec3.norm_sq(v)
  return v.x * v.x + v.y * v.y + v.z * v.z
end

function vec3.norm(v)
  return math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
end

function vec3.unit(v)
  local n = math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
  if n > 1e-15 then
    local inv = 1.0 / n
    return vec3.new(v.x * inv, v.y * inv, v.z * inv)
  end
  return vec3.new(0.0, 0.0, 0.0)
end

function vec3.normalize_inplace(v)
  local n = math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
  if n > 1e-15 then
    local inv = 1.0 / n
    v.x = v.x * inv
    v.y = v.y * inv
    v.z = v.z * inv
  else
    v.x = 0.0
    v.y = 0.0
    v.z = 0.0
  end
  return v
end

function vec3.dist_sq(a, b)
  local dx = a.x - b.x
  local dy = a.y - b.y
  local dz = a.z - b.z
  return dx * dx + dy * dy + dz * dz
end

function vec3.dist(a, b)
  return math.sqrt(vec3.dist_sq(a, b))
end

function vec3.angle(a, b)
  local na = vec3.norm(a)
  local nb = vec3.norm(b)
  if na < 1e-15 or nb < 1e-15 then return 0.0 end
  local c = vec3.dot(a, b) / (na * nb)
  if c > 1.0 then c = 1.0 elseif c < -1.0 then c = -1.0 end
  return math.acos(c)
end

function vec3.project(v, onto)
  local d = vec3.dot(onto, onto)
  if d < 1e-15 then return vec3.zero() end
  local s = vec3.dot(v, onto) / d
  return vec3.new(onto.x * s, onto.y * s, onto.z * s)
end

function vec3.reject(v, from)
  local p = vec3.project(v, from)
  return vec3.new(v.x - p.x, v.y - p.y, v.z - p.z)
end

function vec3.reflect(v, normal)
  local d = 2.0 * vec3.dot(v, normal)
  return vec3.new(v.x - d * normal.x, v.y - d * normal.y, v.z - d * normal.z)
end

function vec3.lerp(a, b, t)
  return vec3.new(
    a.x + (b.x - a.x) * t,
    a.y + (b.y - a.y) * t,
    a.z + (b.z - a.z) * t
  )
end

function vec3.slerp(a, b, t)
  local na = vec3.norm(a)
  local nb = vec3.norm(b)
  if na < 1e-15 or nb < 1e-15 then
    return vec3.lerp(a, b, t)
  end
  local u1 = vec3.new(a.x / na, a.y / na, a.z / na)
  local u2 = vec3.new(b.x / nb, b.y / nb, b.z / nb)
  local dot = vec3.dot(u1, u2)
  if dot > 0.9995 then
    local res = vec3.lerp(a, b, t)
    return res
  end
  if dot < -1.0 then dot = -1.0 elseif dot > 1.0 then dot = 1.0 end
  local theta = math.acos(dot)
  local sin_theta = math.sin(theta)
  local s1 = math.sin((1.0 - t) * theta) / sin_theta
  local s2 = math.sin(t * theta) / sin_theta
  local len = na + (nb - na) * t
  return vec3.new(
    (u1.x * s1 + u2.x * s2) * len,
    (u1.y * s1 + u2.y * s2) * len,
    (u1.z * s1 + u2.z * s2) * len
  )
end

function vec3.skew(v)
  -- Returns 3x3 skew-symmetric matrix (hat map) as array of 3 rows
  return {
    { 0.0, -v.z,  v.y },
    {  v.z, 0.0, -v.x },
    { -v.y,  v.x, 0.0 }
  }
end

-- Metamethods
function vec3_mt.__add(a, b)
  return vec3.new(a.x + b.x, a.y + b.y, a.z + b.z)
end

function vec3_mt.__sub(a, b)
  return vec3.new(a.x - b.x, a.y - b.y, a.z - b.z)
end

function vec3_mt.__mul(a, b)
  if type(a) == "number" then
    return vec3.new(b.x * a, b.y * a, b.z * a)
  elseif type(b) == "number" then
    return vec3.new(a.x * b, a.y * b, a.z * b)
  else
    -- Component-wise multiplication
    return vec3.new(a.x * b.x, a.y * b.y, a.z * b.z)
  end
end

function vec3_mt.__div(a, b)
  if type(b) == "number" then
    local inv = 1.0 / b
    return vec3.new(a.x * inv, a.y * inv, a.z * inv)
  else
    return vec3.new(a.x / b.x, a.y / b.y, a.z / b.z)
  end
end

function vec3_mt.__unm(v)
  return vec3.new(-v.x, -v.y, -v.z)
end

function vec3_mt.__eq(a, b)
  return math.abs(a.x - b.x) < 1e-12 and
         math.abs(a.y - b.y) < 1e-12 and
         math.abs(a.z - b.z) < 1e-12
end

function vec3_mt.__tostring(v)
  return string.format("vec3(%.6f, %.6f, %.6f)", v.x, v.y, v.z)
end

return vec3
