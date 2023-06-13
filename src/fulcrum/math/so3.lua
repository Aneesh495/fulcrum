local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")

local so3 = {}

function so3.hat(v)
  -- Maps vector v in R^3 to skew-symmetric matrix in so(3)
  return mat3.new(
     0.0, -v.z,  v.y,
     v.z,  0.0, -v.x,
    -v.y,  v.x,  0.0
  )
end

function so3.vee(S)
  -- Maps skew-symmetric matrix in so(3) to vector in R^3
  return vec3.new(S[3][2], S[1][3], S[2][1])
end

function so3.bracket(A, B)
  -- Lie bracket [A, B] = A*B - B*A
  return (A * B) - (B * A)
end

function so3.exp(phi)
  -- Matrix exponential from so(3) vector to SO(3) rotation matrix
  -- Rodrigues formula with high-order series expansions near zero
  local theta_sq = vec3.norm_sq(phi)
  local theta = math.sqrt(theta_sq)
  local K = so3.hat(phi)
  local K2 = K * K

  local a, b
  if theta_sq < 1e-8 then
    -- Taylor expansion near theta = 0
    a = 1.0 - theta_sq / 6.0 + (theta_sq * theta_sq) / 120.0
    b = 0.5 - theta_sq / 24.0 + (theta_sq * theta_sq) / 720.0
  else
    a = math.sin(theta) / theta
    b = (1.0 - math.cos(theta)) / theta_sq
  end

  local I = mat3.identity()
  return I + (K * a) + (K2 * b)
end

function so3.log(R)
  -- Matrix logarithm from SO(3) to vector in R^3 (so(3) algebra)
  local tr = mat3.trace(R)
  local cos_theta = math.max(-1.0, math.min(1.0, (tr - 1.0) * 0.5))
  local theta = math.acos(cos_theta)

  if theta < 1e-7 then
    -- Near identity: Taylor expansion
    local a = 0.5 + (theta * theta) / 12.0
    return vec3.new(
      (R[3][2] - R[2][3]) * a,
      (R[1][3] - R[3][1]) * a,
      (R[2][1] - R[1][2]) * a
    )
  elseif math.abs(theta - math.pi) < 1e-4 then
    -- Near pi singularity: use diagonal components
    local xx = (R[1][1] + 1.0) * 0.5
    local yy = (R[2][2] + 1.0) * 0.5
    local zz = (R[3][3] + 1.0) * 0.5
    local xy = (R[1][2] + R[2][1]) * 0.25
    local xz = (R[1][3] + R[3][1]) * 0.25
    local yz = (R[2][3] + R[3][2]) * 0.25

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
    local unit_axis = vec3.unit(vec3.new(rx, ry, rz))
    return unit_axis * theta
  else
    local factor = theta / (2.0 * math.sin(theta))
    return vec3.new(
      (R[3][2] - R[2][3]) * factor,
      (R[1][3] - R[3][1]) * factor,
      (R[2][1] - R[1][2]) * factor
    )
  end
end

function so3.left_jacobian(phi)
  -- J_l(phi) = (sin theta / theta)*I + (1 - sin theta / theta)*(phi*phi^T / theta^2) + ((1 - cos theta)/theta^2)*hat(phi)
  local theta_sq = vec3.norm_sq(phi)
  local theta = math.sqrt(theta_sq)
  local I = mat3.identity()
  local K = so3.hat(phi)

  if theta_sq < 1e-8 then
    return I + (K * 0.5)
  end

  local sin_t = math.sin(theta)
  local cos_t = math.cos(theta)
  local a = sin_t / theta
  local b = (1.0 - cos_t) / theta_sq
  local c = (1.0 - a) / theta_sq

  local phi_phi_T = mat3.new(
    phi.x * phi.x, phi.x * phi.y, phi.x * phi.z,
    phi.y * phi.x, phi.y * phi.y, phi.y * phi.z,
    phi.z * phi.x, phi.z * phi.y, phi.z * phi.z
  )

  return (I * a) + (K * b) + (phi_phi_T * c)
end

function so3.right_jacobian(phi)
  -- J_r(phi) = J_l(-phi)
  return so3.left_jacobian(-phi)
end

function so3.left_jacobian_inv(phi)
  -- Inverse of left Jacobian J_l^{-1}(phi)
  -- = (theta / 2 * cot(theta / 2))*I + (1 - theta / 2 * cot(theta / 2))*(phi*phi^T / theta^2) - 0.5 * hat(phi)
  local theta_sq = vec3.norm_sq(phi)
  local theta = math.sqrt(theta_sq)
  local I = mat3.identity()
  local K = so3.hat(phi)

  local alpha
  if theta_sq < 1e-8 then
    alpha = 1.0 - theta_sq / 12.0 - (theta_sq * theta_sq) / 720.0
  else
    local half_theta = 0.5 * theta
    alpha = half_theta / math.tan(half_theta)
  end

  local phi_phi_T = mat3.new(
    phi.x * phi.x, phi.x * phi.y, phi.x * phi.z,
    phi.y * phi.x, phi.y * phi.y, phi.y * phi.z,
    phi.z * phi.x, phi.z * phi.y, phi.z * phi.z
  )
  local factor = (1.0 - alpha) / (theta_sq < 1e-12 and 1.0 or theta_sq)

  return (I * alpha) - (K * 0.5) + (phi_phi_T * factor)
end

function so3.right_jacobian_inv(phi)
  -- J_r^{-1}(phi) = J_l^{-1}(-phi)
  return so3.left_jacobian_inv(-phi)
end

function so3.adjoint(R)
  -- Adjoint representation for SO(3) is the matrix R itself
  return mat3.clone(R)
end

function so3.distance(R1, R2)
  -- Geodesic distance on SO(3) manifold: ||log(R1^T * R2)||
  local R1t = mat3.transpose(R1)
  local diff = R1t * R2
  return vec3.norm(so3.log(diff))
end

return so3
