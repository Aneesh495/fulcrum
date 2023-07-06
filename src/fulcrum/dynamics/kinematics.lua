local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")

local kinematics = {}

function kinematics.forward_kinematics(m, q)
  -- Computes:
  -- X_tree[i]: spatial transform from parent of body i to body i: {}^{p(i)} X_i
  -- X_world[i]: spatial transform from world (frame 0) to body i: {}^0 X_i
  local X_tree = {}
  local X_world = {}

  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local qi = {}
    if q_idx and q_idx.dof > 0 then
      for k = q_idx.start, q_idx.finish do
        table.insert(qi, q[k] or 0.0)
      end
    end

    local X_J = b.joint.transform(qi)
    -- Transform across joint from parent to child: X_J * X_T
    local X_pi = X_J * b.X_T
    X_tree[i] = X_pi

    if b.parent_id == 0 then
      X_world[i] = X_pi
    else
      -- X_0i = X_pi * X_0p
      X_world[i] = X_pi * X_world[b.parent_id]
    end
  end

  return X_tree, X_world
end

function kinematics.forward_velocities(m, q, q_dot)
  local X_tree, X_world = kinematics.forward_kinematics(m, q)
  local v = {}

  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local qi = {}
    local qd_i = {}

    if q_idx and q_idx.dof > 0 then
      for k = q_idx.start, q_idx.finish do
        table.insert(qi, q[k] or 0.0)
        table.insert(qd_i, q_dot[k] or 0.0)
      end
    end

    local S = b.joint.subspace(qi)
    local v_joint = matrix.mul_vec(S, qd_i)

    if b.parent_id == 0 then
      v[i] = v_joint
    else
      -- v_i = X_pi * v_p + v_joint
      local v_from_parent = spatial.transform_motion(X_tree[i], v[b.parent_id])
      v[i] = spatial.add(v_from_parent, v_joint)
    end
  end

  return v, X_tree, X_world
end

function kinematics.center_of_mass(m, q)
  local _, X_world = kinematics.forward_kinematics(m, q)
  local total_mass = 0.0
  local com_numerator = vec3.zero()

  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    total_mass = total_mass + b.mass

    -- Extract body world position and rotation from X_world
    -- X_0i = [ R, 0 ; -R*skew(p), R ]
    -- Note: X_world maps from world to body frame: v_body = X * v_world
    -- The origin of body frame in world is p_i = -R^T * (-R*skew(p))_extracted
    -- More directly: invert X_world to get transform from body to world
    local X_inv = matrix.inv(X_world[i])
    -- Origin of body i in world frame:
    local R_wb = mat3.new(
      X_inv[1][1], X_inv[1][2], X_inv[1][3],
      X_inv[2][1], X_inv[2][2], X_inv[2][3],
      X_inv[3][1], X_inv[3][2], X_inv[3][3]
    )
    -- Translational component of X_inv is at bottom-left: -R_wb * skew(p_body)
    -- For point in body frame: p_w = R_wb * b.com + p_body_origin
    -- Let's extract position by transforming [0,0,0]:
    -- An easy way: spatial motion of pure linear displacement
    local px = -X_inv[5][3]
    local py = X_inv[4][3]
    local pz = -X_inv[4][2]
    local p_origin = vec3.new(px, py, pz)
    local body_com_world = p_origin + (R_wb * b.com)

    com_numerator = com_numerator + (body_com_world * b.mass)
  end

  local system_com = com_numerator / (total_mass > 0 and total_mass or 1.0)
  return system_com, total_mass
end

return kinematics
