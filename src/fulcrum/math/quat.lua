local vec3 = require("fulcrum.math.vec3")

local quat = {}
local quat_mt = { __index = quat }

function quat.new(w, x, y, z)
  local q = {
    w = w or 1.0,
    x = x or 0.0,
    y = y or 0.0,
    z = z or 0.0
  }
  return setmetatable(q, quat_mt)
end

function quat.identity()
  return quat.new(1.0, 0.0, 0.0, 0.0)
end

function quat.clone(q)
  return quat.new(q.w, q.x, q.y, q.z)
end

function quat.norm_sq(q)
  return q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z
end

function quat.norm(q)
  return math.sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z)
end

function quat.unit(q)
  local n = math.sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z)
  if n > 1e-15 then
    local inv = 1.0 / n
    return quat.new(q.w * inv, q.x * inv, q.y * inv, q.z * inv)
  end
  return quat.identity()
end

function quat.normalize_inplace(q)
  local n = math.sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z)
  if n > 1e-15 then
    local inv = 1.0 / n
    q.w = q.w * inv
    q.x = q.x * inv
    q.y = q.y * inv
    q.z = q.z * inv
  else
    q.w = 1.0
    q.x = 0.0
    q.y = 0.0
    q.z = 0.0
  end
  return q
end

function quat.conjugate(q)
  return quat.new(q.w, -q.x, -q.y, -q.z)
end

function quat.inverse(q)
  local nsq = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z
  if nsq > 1e-15 then
    local inv = 1.0 / nsq
    return quat.new(q.w * inv, -q.x * inv, -q.y * inv, -q.z * inv)
  end
  return quat.identity()
end

function quat.from_axis_angle(axis, angle)
  local na = math.sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z)
  if na < 1e-15 then
    return quat.identity()
  end
  local half = 0.5 * angle
  local s = math.sin(half) / na
  return quat.new(math.cos(half), axis.x * s, axis.y * s, axis.z * s)
end

function quat.to_axis_angle(q)
  local qn = quat.unit(q)
  if qn.w < 0 then
    qn = quat.new(-qn.w, -qn.x, -qn.y, -qn.z)
  end
  local angle = 2.0 * math.acos(math.max(-1.0, math.min(1.0, qn.w)))
  local s = math.sqrt(1.0 - qn.w * qn.w)
  local axis
  if s < 1e-7 then
    axis = vec3.new(1.0, 0.0, 0.0)
  else
    local inv = 1.0 / s
    axis = vec3.new(qn.x * inv, qn.y * inv, qn.z * inv)
  end
  return axis, angle
end

function quat.from_euler(roll, pitch, yaw)
  local cr = math.cos(roll * 0.5)
  local sr = math.sin(roll * 0.5)
  local cp = math.cos(pitch * 0.5)
  local sp = math.sin(pitch * 0.5)
  local cy = math.cos(yaw * 0.5)
  local sy = math.sin(yaw * 0.5)

  return quat.new(
    cr * cp * cy + sr * sp * sy,
    sr * cp * cy - cr * sp * sy,
    cr * sp * cy + sr * cp * sy,
    cr * cp * sy - sr * sp * cy
  )
end

function quat.to_euler(q)
  local qn = quat.unit(q)
  -- Roll (x-axis rotation)
  local sinr_cosp = 2.0 * (qn.w * qn.x + qn.y * qn.z)
  local cosr_cosp = 1.0 - 2.0 * (qn.x * qn.x + qn.y * qn.y)
  local roll = math.atan(sinr_cosp, cosr_cosp)

  -- Pitch (y-axis rotation)
  local sinp = 2.0 * (qn.w * qn.y - qn.z * qn.x)
  local pitch
  if math.abs(sinp) >= 1.0 then
    pitch = (sinp >= 0) and (math.pi * 0.5) or (-math.pi * 0.5)
  else
    pitch = math.asin(sinp)
  end

  -- Yaw (z-axis rotation)
  local siny_cosp = 2.0 * (qn.w * qn.z + qn.x * qn.y)
  local cosy_cosp = 1.0 - 2.0 * (qn.y * qn.y + qn.z * qn.z)
  local yaw = math.atan(siny_cosp, cosy_cosp)

  return roll, pitch, yaw
end

function quat.from_rotation_matrix(R)
  -- Shepperd algorithm for robust quaternion extraction
  local tr = R[1][1] + R[2][2] + R[3][3]
  local q = quat.new()

  if tr > 0.0 then
    local s = 2.0 * math.sqrt(tr + 1.0)
    local inv = 1.0 / s
    q.w = 0.25 * s
    q.x = (R[3][2] - R[2][3]) * inv
    q.y = (R[1][3] - R[3][1]) * inv
    q.z = (R[2][1] - R[1][2]) * inv
  elseif (R[1][1] > R[2][2]) and (R[1][1] > R[3][3]) then
    local s = 2.0 * math.sqrt(1.0 + R[1][1] - R[2][2] - R[3][3])
    local inv = 1.0 / s
    q.w = (R[3][2] - R[2][3]) * inv
    q.x = 0.25 * s
    q.y = (R[1][2] + R[2][1]) * inv
    q.z = (R[1][3] + R[3][1]) * inv
  elseif R[2][2] > R[3][3] then
    local s = 2.0 * math.sqrt(1.0 + R[2][2] - R[1][1] - R[3][3])
    local inv = 1.0 / s
    q.w = (R[1][3] - R[3][1]) * inv
    q.x = (R[1][2] + R[2][1]) * inv
    q.y = 0.25 * s
    q.z = (R[2][3] + R[3][2]) * inv
  else
    local s = 2.0 * math.sqrt(1.0 + R[3][3] - R[1][1] - R[2][2])
    local inv = 1.0 / s
    q.w = (R[2][1] - R[1][2]) * inv
    q.x = (R[1][3] + R[3][1]) * inv
    q.y = (R[2][3] + R[3][2]) * inv
    q.z = 0.25 * s
  end

  return quat.unit(q)
end

function quat.to_rotation_matrix(q)
  local qn = quat.unit(q)
  local w, x, y, z = qn.w, qn.x, qn.y, qn.z

  local xx = x * x
  local yy = y * y
  local zz = z * z
  local xy = x * y
  local xz = x * z
  local yz = y * z
  local wx = w * x
  local wy = w * y
  local wz = w * z

  return {
    { 1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy) },
    { 2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx) },
    { 2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy) }
  }
end

function quat.rotate_vector(q, v)
  -- Vector rotation: v' = v + 2*w*(qv x v) + 2*(qv x (qv x v))
  local qv = vec3.new(q.x, q.y, q.z)
  local t = vec3.cross(qv, v) * 2.0
  local cross_t = vec3.cross(qv, t)
  return v + (t * q.w) + cross_t
end

function quat.slerp(q1, q2, t)
  local dot = q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z
  local w2, x2, y2, z2 = q2.w, q2.x, q2.y, q2.z

  -- Antipodal flip for shortest path
  if dot < 0.0 then
    dot = -dot
    w2, x2, y2, z2 = -w2, -x2, -y2, -z2
  end

  if dot > 0.9995 then
    -- Linear interpolation for very close orientations
    local q_res = quat.new(
      q1.w + (w2 - q1.w) * t,
      q1.x + (x2 - q1.x) * t,
      q1.y + (y2 - q1.y) * t,
      q1.z + (z2 - q1.z) * t
    )
    return quat.unit(q_res)
  end

  dot = math.max(-1.0, math.min(1.0, dot))
  local theta_0 = math.acos(dot)
  local theta = theta_0 * t
  local sin_theta = math.sin(theta)
  local sin_theta_0 = math.sin(theta_0)

  local s1 = math.cos(theta) - dot * sin_theta / sin_theta_0
  local s2 = sin_theta / sin_theta_0

  return quat.new(
    q1.w * s1 + w2 * s2,
    q1.x * s1 + x2 * s2,
    q1.y * s1 + y2 * s2,
    q1.z * s1 + z2 * s2
  )
end

function quat.derivative(q, omega, is_body_fixed)
  -- Time derivative:
  -- If body-fixed: q_dot = 0.5 * q * [0, omega]
  -- If spatial:    q_dot = 0.5 * [0, omega] * q
  local wq = quat.new(0.0, omega.x, omega.y, omega.z)
  if is_body_fixed ~= false then
    return (q * wq) * 0.5
  else
    return (wq * q) * 0.5
  end
end

function quat.integrate_angular_velocity(q, omega, dt)
  local norm_omega = vec3.norm(omega)
  if norm_omega < 1e-12 then
    return quat.clone(q)
  end
  local half_angle = 0.5 * norm_omega * dt
  local s = math.sin(half_angle) / norm_omega
  local delta_q = quat.new(
    math.cos(half_angle),
    omega.x * s,
    omega.y * s,
    omega.z * s
  )
  return quat.unit(q * delta_q)
end

function quat.exp(v)
  local theta = vec3.norm(v)
  if theta < 1e-12 then
    return quat.identity()
  end
  local s = math.sin(theta) / theta
  return quat.new(math.cos(theta), v.x * s, v.y * s, v.z * s)
end

function quat.log(q)
  local qn = quat.unit(q)
  local v_norm = math.sqrt(qn.x * qn.x + qn.y * qn.y + qn.z * qn.z)
  if v_norm < 1e-12 then
    return vec3.zero()
  end
  local theta = math.atan(v_norm, qn.w)
  local s = theta / v_norm
  return vec3.new(qn.x * s, qn.y * s, qn.z * s)
end

-- Metamethods
function quat_mt.__mul(a, b)
  if type(a) == "number" then
    return quat.new(b.w * a, b.x * a, b.y * a, b.z * a)
  elseif type(b) == "number" then
    return quat.new(a.w * b, a.x * b, a.y * b, a.z * b)
  else
    -- Hamilton product
    return quat.new(
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
      a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
      a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
      a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    )
  end
end

function quat_mt.__add(a, b)
  return quat.new(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function quat_mt.__sub(a, b)
  return quat.new(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

function quat_mt.__unm(q)
  return quat.new(-q.w, -q.x, -q.y, -q.z)
end

function quat_mt.__eq(a, b)
  local dot = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z
  return math.abs(math.abs(dot) - 1.0) < 1e-10
end

function quat_mt.__tostring(q)
  return string.format("quat(%.6f; %.6f, %.6f, %.6f)", q.w, q.x, q.y, q.z)
end

return quat
