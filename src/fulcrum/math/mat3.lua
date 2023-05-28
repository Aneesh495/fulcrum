local vec3 = require("fulcrum.math.vec3")

local mat3 = {}
local mat3_mt = { __index = mat3 }

function mat3.new(m11, m12, m13, m21, m22, m23, m31, m32, m33)
  local m = {
    { m11 or 0.0, m12 or 0.0, m13 or 0.0 },
    { m21 or 0.0, m22 or 0.0, m23 or 0.0 },
    { m31 or 0.0, m32 or 0.0, m33 or 0.0 }
  }
  return setmetatable(m, mat3_mt)
end

function mat3.identity()
  return mat3.new(
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  )
end

function mat3.zero()
  return mat3.new(
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0
  )
end

function mat3.from_diag(dx, dy, dz)
  return mat3.new(
    dx or 0.0, 0.0, 0.0,
    0.0, dy or 0.0, 0.0,
    0.0, 0.0, dz or 0.0
  )
end

function mat3.from_rows(r1, r2, r3)
  return mat3.new(
    r1.x or r1[1], r1.y or r1[2], r1.z or r1[3],
    r2.x or r2[1], r2.y or r2[2], r2.z or r2[3],
    r3.x or r3[1], r3.y or r3[2], r3.z or r3[3]
  )
end

function mat3.from_cols(c1, c2, c3)
  return mat3.new(
    c1.x or c1[1], c2.x or c2[1], c3.x or c3[1],
    c1.y or c1[2], c2.y or c2[2], c3.y or c3[2],
    c1.z or c1[3], c2.z or c2[3], c3.z or c3[3]
  )
end

function mat3.clone(m)
  return mat3.new(
    m[1][1], m[1][2], m[1][3],
    m[2][1], m[2][2], m[2][3],
    m[3][1], m[3][2], m[3][3]
  )
end

function mat3.transpose(m)
  return mat3.new(
    m[1][1], m[2][1], m[3][1],
    m[1][2], m[2][2], m[3][2],
    m[1][3], m[2][3], m[3][3]
  )
end

function mat3.trace(m)
  return m[1][1] + m[2][2] + m[3][3]
end

function mat3.det(m)
  local m11, m12, m13 = m[1][1], m[1][2], m[1][3]
  local m21, m22, m23 = m[2][1], m[2][2], m[2][3]
  local m31, m32, m33 = m[3][1], m[3][2], m[3][3]

  return m11 * (m22 * m33 - m23 * m32)
       - m12 * (m21 * m33 - m23 * m31)
       + m13 * (m21 * m32 - m22 * m31)
end

function mat3.inv(m)
  local m11, m12, m13 = m[1][1], m[1][2], m[1][3]
  local m21, m22, m23 = m[2][1], m[2][2], m[2][3]
  local m31, m32, m33 = m[3][1], m[3][2], m[3][3]

  local c11 =  (m22 * m33 - m23 * m32)
  local c12 = -(m21 * m33 - m23 * m31)
  local c13 =  (m21 * m32 - m22 * m31)

  local det = m11 * c11 + m12 * c12 + m13 * c13
  if math.abs(det) < 1e-15 then
    error("Matrix is singular and cannot be inverted")
  end

  local inv_det = 1.0 / det
  local c21 = -(m12 * m33 - m13 * m32)
  local c22 =  (m11 * m33 - m13 * m31)
  local c23 = -(m11 * m32 - m12 * m31)

  local c31 =  (m12 * m23 - m13 * m22)
  local c32 = -(m11 * m23 - m13 * m21)
  local c33 =  (m11 * m22 - m12 * m21)

  return mat3.new(
    c11 * inv_det, c21 * inv_det, c31 * inv_det,
    c12 * inv_det, c22 * inv_det, c32 * inv_det,
    c13 * inv_det, c23 * inv_det, c33 * inv_det
  )
end

function mat3.mul(a, b)
  local res = mat3.zero()
  for i = 1, 3 do
    for j = 1, 3 do
      res[i][j] = a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j]
    end
  end
  return res
end

function mat3.mul_vec(m, v)
  return vec3.new(
    m[1][1] * v.x + m[1][2] * v.y + m[1][3] * v.z,
    m[2][1] * v.x + m[2][2] * v.y + m[2][3] * v.z,
    m[3][1] * v.x + m[3][2] * v.y + m[3][3] * v.z
  )
end

function mat3.from_axis_angle(axis, angle)
  local na = math.sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z)
  if na < 1e-15 then
    return mat3.identity()
  end
  local x = axis.x / na
  local y = axis.y / na
  local z = axis.z / na

  local c = math.cos(angle)
  local s = math.sin(angle)
  local c1 = 1.0 - c

  return mat3.new(
    c + x * x * c1,     x * y * c1 - z * s, x * z * c1 + y * s,
    y * x * c1 + z * s, c + y * y * c1,     y * z * c1 - x * s,
    z * x * c1 - y * s, z * y * c1 + x * s, c + z * z * c1
  )
end

function mat3.to_axis_angle(m)
  local tr = mat3.trace(m)
  local cos_theta = math.max(-1.0, math.min(1.0, (tr - 1.0) * 0.5))
  local angle = math.acos(cos_theta)

  if angle < 1e-7 then
    return vec3.new(1.0, 0.0, 0.0), 0.0
  elseif math.abs(angle - math.pi) < 1e-4 then
    -- Near 180 degrees singularity
    local xx = (m[1][1] + 1.0) * 0.5
    local yy = (m[2][2] + 1.0) * 0.5
    local zz = (m[3][3] + 1.0) * 0.5
    local xy = (m[1][2] + m[2][1]) * 0.25
    local xz = (m[1][3] + m[3][1]) * 0.25
    local yz = (m[2][3] + m[3][2]) * 0.25

    local rx = math.sqrt(math.max(0.0, xx))
    local ry = math.sqrt(math.max(0.0, yy))
    local rz = math.sqrt(math.max(0.0, zz))

    if rx >= ry and rx >= rz then
      if xy < 0 then ry = -ry end
      if xz < 0 then rz = -rz end
    elseif ry >= rz then
      if xy < 0 then rx = -rx end
      if yz < 0 then rz = -rz end
    else
      if xz < 0 then rx = -rx end
      if yz < 0 then ry = -ry end
    end
    return vec3.unit(vec3.new(rx, ry, rz)), angle
  else
    local s = 2.0 * math.sin(angle)
    local inv = 1.0 / s
    local axis = vec3.new(
      (m[3][2] - m[2][3]) * inv,
      (m[1][3] - m[3][1]) * inv,
      (m[2][1] - m[1][2]) * inv
    )
    return vec3.unit(axis), angle
  end
end

function mat3.orthonormalize(m)
  -- Gram-Schmidt process on columns
  local c1 = vec3.new(m[1][1], m[2][1], m[3][1])
  local c2 = vec3.new(m[1][2], m[2][2], m[3][2])
  local c3 = vec3.new(m[1][3], m[2][3], m[3][3])

  local u1 = vec3.unit(c1)
  local u2 = vec3.unit(c2 - vec3.project(c2, u1))
  local u3 = vec3.cross(u1, u2)

  return mat3.from_cols(u1, u2, u3)
end

-- Metamethods
function mat3_mt.__add(a, b)
  local res = mat3.zero()
  for i = 1, 3 do
    for j = 1, 3 do
      res[i][j] = a[i][j] + b[i][j]
    end
  end
  return res
end

function mat3_mt.__sub(a, b)
  local res = mat3.zero()
  for i = 1, 3 do
    for j = 1, 3 do
      res[i][j] = a[i][j] - b[i][j]
    end
  end
  return res
end

function mat3_mt.__mul(a, b)
  if type(a) == "number" then
    local res = mat3.zero()
    for i = 1, 3 do
      for j = 1, 3 do
        res[i][j] = b[i][j] * a
      end
    end
    return res
  elseif type(b) == "number" then
    local res = mat3.zero()
    for i = 1, 3 do
      for j = 1, 3 do
        res[i][j] = a[i][j] * b
      end
    end
    return res
  elseif b.x and b.y and b.z then
    -- Matrix times vector
    return mat3.mul_vec(a, b)
  else
    -- Matrix times matrix
    return mat3.mul(a, b)
  end
end

function mat3_mt.__unm(m)
  local res = mat3.zero()
  for i = 1, 3 do
    for j = 1, 3 do
      res[i][j] = -m[i][j]
    end
  end
  return res
end

function mat3_mt.__eq(a, b)
  for i = 1, 3 do
    for j = 1, 3 do
      if math.abs(a[i][j] - b[i][j]) > 1e-10 then
        return false
      end
    end
  end
  return true
end

function mat3_mt.__tostring(m)
  return string.format("mat3(\n  [%.5f, %.5f, %.5f]\n  [%.5f, %.5f, %.5f]\n  [%.5f, %.5f, %.5f]\n)",
    m[1][1], m[1][2], m[1][3],
    m[2][1], m[2][2], m[2][3],
    m[3][1], m[3][2], m[3][3]
  )
end

return mat3
