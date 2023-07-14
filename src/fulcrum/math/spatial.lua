local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")

local spatial = {}

function spatial.vector(w, v)
  -- 6D spatial motion or force vector: elements 1..3 angular, 4..6 linear
  return {
    w.x or w[1] or 0.0,
    w.y or w[2] or 0.0,
    w.z or w[3] or 0.0,
    v.x or v[1] or 0.0,
    v.y or v[2] or 0.0,
    v.z or v[3] or 0.0
  }
end

function spatial.zero()
  return { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
end

function spatial.angular(s)
  return vec3.new(s[1], s[2], s[3])
end

function spatial.linear(s)
  return vec3.new(s[4], s[5], s[6])
end

function spatial.add(a, b)
  return {
    a[1] + b[1],
    a[2] + b[2],
    a[3] + b[3],
    a[4] + b[4],
    a[5] + b[5],
    a[6] + b[6]
  }
end

function spatial.sub(a, b)
  return {
    a[1] - b[1],
    a[2] - b[2],
    a[3] - b[3],
    a[4] - b[4],
    a[5] - b[5],
    a[6] - b[6]
  }
end

function spatial.scale(s, k)
  return {
    s[1] * k,
    s[2] * k,
    s[3] * k,
    s[4] * k,
    s[5] * k,
    s[6] * k
  }
end

function spatial.dot(force, motion)
  -- Scalar product between force vector and motion vector (power)
  return force[1] * motion[1] + force[2] * motion[2] + force[3] * motion[3] +
         force[4] * motion[4] + force[5] * motion[5] + force[6] * motion[6]
end

function spatial.cross_motion(v, m)
  -- Spatial cross product for motion vectors: v x m
  -- w_res = w_v x w_m
  -- v_res = w_v x v_m + v_v x w_m
  local wv = vec3.new(v[1], v[2], v[3])
  local vv = vec3.new(v[4], v[5], v[6])
  local wm = vec3.new(m[1], m[2], m[3])
  local vm = vec3.new(m[4], m[5], m[6])

  local w_res = vec3.cross(wv, wm)
  local v_res = vec3.cross(wv, vm) + vec3.cross(vv, wm)

  return { w_res.x, w_res.y, w_res.z, v_res.x, v_res.y, v_res.z }
end

function spatial.cross_force(v, f)
  -- Spatial cross product for force vectors: v x* f
  -- n_res = w_v x n_f + v_v x f_f
  -- f_res = w_v x f_f
  local wv = vec3.new(v[1], v[2], v[3])
  local vv = vec3.new(v[4], v[5], v[6])
  local nf = vec3.new(f[1], f[2], f[3])
  local ff = vec3.new(f[4], f[5], f[6])

  local n_res = vec3.cross(wv, nf) + vec3.cross(vv, ff)
  local f_res = vec3.cross(wv, ff)

  return { n_res.x, n_res.y, n_res.z, f_res.x, f_res.y, f_res.z }
end

function spatial.transform(R, p)
  -- Spatial coordinate transformation matrix X = [ R, 0 ; -R * skew(p), R ]
  -- Transforms motion vectors from frame B to frame A: v_A = X * v_B
  local X = matrix.zeros(6, 6)
  local p_skew = vec3.skew(p)

  -- Top-left 3x3: R
  for i = 1, 3 do
    for j = 1, 3 do
      X[i][j] = R[i][j]
    end
  end

  -- Bottom-right 3x3: R
  for i = 1, 3 do
    for j = 1, 3 do
      X[3 + i][3 + j] = R[i][j]
    end
  end

  -- Bottom-left 3x3: -R * skew(p)
  for i = 1, 3 do
    for j = 1, 3 do
      local sum = 0.0
      for k = 1, 3 do
        sum = sum + R[i][k] * p_skew[k][j]
      end
      X[3 + i][j] = -sum
    end
  end

  return X
end

function spatial.transform_motion(X, v)
  -- Efficient motion transform without full 6x6 matrix multiply if R and p are provided
  -- Or multiply 6x6 by 6-vector
  return matrix.mul_vec(X, v)
end

function spatial.transform_force(X, f)
  -- Force transform by duality (power invariance):
  -- If v_B = X * v_A, then f_A = X^T * f_B
  local Xt = matrix.transpose(X)
  return matrix.mul_vec(Xt, f)
end

function spatial.inertia(mass, com, I_rot)
  -- 6x6 Spatial rigid body inertia matrix
  -- I = [ I_rot + m * skew(c) * skew(c)^T,  m * skew(c)   ]
  --     [ m * skew(c)^T,                     m * eye(3)    ]
  -- where c is the center of mass vector
  local I_sp = matrix.zeros(6, 6)
  local c_skew = vec3.skew(com)

  -- Bottom-right: m * eye(3)
  for i = 1, 3 do
    I_sp[3 + i][3 + i] = mass
  end

  -- Top-right: m * skew(c)
  for i = 1, 3 do
    for j = 1, 3 do
      I_sp[i][3 + j] = mass * c_skew[i][j]
    end
  end

  -- Bottom-left: -m * skew(c) = m * skew(c)^T
  for i = 1, 3 do
    for j = 1, 3 do
      I_sp[3 + i][j] = -mass * c_skew[i][j]
    end
  end

  -- Top-left: I_rot + m * c~ * c~^T = I_rot - m * c~ * c~
  for i = 1, 3 do
    for j = 1, 3 do
      local cc = 0.0
      for k = 1, 3 do
        cc = cc + c_skew[i][k] * c_skew[k][j]
      end
      I_sp[i][j] = I_rot[i][j] - mass * cc
    end
  end

  return I_sp
end

function spatial.kinetic_energy(I_sp, v)
  -- T = 0.5 * v^T * I * v
  local Iv = matrix.mul_vec(I_sp, v)
  local e = 0.0
  for i = 1, 6 do
    e = e + v[i] * Iv[i]
  end
  return 0.5 * e
end

return spatial
