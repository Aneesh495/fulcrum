local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local so3 = require("fulcrum.math.so3")

local se3 = {}
local se3_mt = { __index = se3 }

function se3.new(R, p)
  local t = {
    R = R or mat3.identity(),
    p = p or vec3.zero()
  }
  return setmetatable(t, se3_mt)
end

function se3.identity()
  return se3.new(mat3.identity(), vec3.zero())
end

function se3.from_matrix4(M)
  local R = mat3.new(
    M[1][1], M[1][2], M[1][3],
    M[2][1], M[2][2], M[2][3],
    M[3][1], M[3][2], M[3][3]
  )
  local p = vec3.new(M[1][4], M[2][4], M[3][4])
  return se3.new(R, p)
end

function se3.to_matrix4(T)
  return {
    { T.R[1][1], T.R[1][2], T.R[1][3], T.p.x },
    { T.R[2][1], T.R[2][2], T.R[2][3], T.p.y },
    { T.R[3][1], T.R[3][2], T.R[3][3], T.p.z },
    { 0.0,       0.0,       0.0,       1.0 }
  }
end

function se3.inv(T)
  local Rt = mat3.transpose(T.R)
  local inv_p = -(Rt * T.p)
  return se3.new(Rt, inv_p)
end

function se3.mul(A, B)
  -- T_AB * T_BC = [ R_A * R_B,  p_A + R_A * p_B ]
  local R = A.R * B.R
  local p = A.p + (A.R * B.p)
  return se3.new(R, p)
end

function se3.act(T, v)
  -- Transform a 3D point: p' = R * v + p
  return (T.R * v) + T.p
end

function se3.act_vector(T, v)
  -- Transform a 3D direction vector: v' = R * v
  return T.R * v
end

function se3.hat(twist)
  -- Twist xi = (v, w): elements 1..3 linear, 4..6 angular
  local v = vec3.new(twist[1], twist[2], twist[3])
  local w = vec3.new(twist[4], twist[5], twist[6])
  local w_hat = so3.hat(w)

  return {
    { w_hat[1][1], w_hat[1][2], w_hat[1][3], v.x },
    { w_hat[2][1], w_hat[2][2], w_hat[2][3], v.y },
    { w_hat[3][1], w_hat[3][2], w_hat[3][3], v.z },
    { 0.0,         0.0,         0.0,         0.0 }
  }
end

function se3.vee(Xi)
  local w = so3.vee(Xi)
  local v = vec3.new(Xi[1][4], Xi[2][4], Xi[3][4])
  return { v.x, v.y, v.z, w.x, w.y, w.z }
end

function se3.exp(twist)
  -- Exponential map: SE(3) = exp(hat(twist))
  -- twist = (v, w)
  local v = vec3.new(twist[1], twist[2], twist[3])
  local w = vec3.new(twist[4], twist[5], twist[6])

  local theta_sq = vec3.norm_sq(w)
  local R = so3.exp(w)
  local p

  if theta_sq < 1e-12 then
    p = v
  else
    local Jl = so3.left_jacobian(w)
    p = Jl * v
  end

  return se3.new(R, p)
end

function se3.log(T)
  -- Logarithmic map: twist = log(T) in R^6
  local w = so3.log(T.R)
  local theta_sq = vec3.norm_sq(w)
  local v

  if theta_sq < 1e-12 then
    v = T.p
  else
    local Jl_inv = so3.left_jacobian_inv(w)
    v = Jl_inv * T.p
  end

  return { v.x, v.y, v.z, w.x, w.y, w.z }
end

function se3.adjoint(T)
  -- 6x6 Adjoint matrix Ad_T
  -- For twist (v, w):
  -- Ad_T = [ R,  hat(p) * R ]
  --        [ 0,  R          ]
  local Ad = matrix.zeros(6, 6)
  local p_skew = so3.hat(T.p)
  local p_skew_R = p_skew * T.R

  for i = 1, 3 do
    for j = 1, 3 do
      Ad[i][j] = T.R[i][j]
      Ad[i][3 + j] = p_skew_R[i][j]
      Ad[3 + i][3 + j] = T.R[i][j]
    end
  end

  return Ad
end

function se3.interpolate(T1, T2, t)
  -- Geodesic interpolation along screw axis
  local diff = se3.inv(T1) * T2
  local twist = se3.log(diff)
  local scaled_twist = {
    twist[1] * t,
    twist[2] * t,
    twist[3] * t,
    twist[4] * t,
    twist[5] * t,
    twist[6] * t
  }
  return T1 * se3.exp(scaled_twist)
end

-- Metamethods
function se3_mt.__mul(A, B)
  return se3.mul(A, B)
end

function se3_mt.__tostring(T)
  return string.format("se3(p=[%.4f, %.4f, %.4f])", T.p.x, T.p.y, T.p.z)
end

return se3
