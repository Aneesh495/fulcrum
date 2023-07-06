local vec3 = require("fulcrum.math.vec3")
local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")
local spatial = require("fulcrum.math.spatial")
local kinematics = require("fulcrum.dynamics.kinematics")

local jacobians = {}

function jacobians.geometric_jacobian(m, q, body_id, in_world_frame)
  -- Computes 6 x dof geometric Jacobian for body_id
  -- If in_world_frame is true, expresses spatial velocity in world coordinates
  local X_tree, X_world = kinematics.forward_kinematics(m, q)
  local J = matrix.zeros(6, m.dof)

  -- Trace ancestors of body_id up to root
  local curr = body_id
  local ancestors = {}
  while curr > 0 do
    ancestors[curr] = true
    curr = m.bodies[curr].parent_id
  end

  for j = 1, m.num_bodies do
    if ancestors[j] then
      local b = m.bodies[j]
      local q_idx = m.q_indices[j]
      local qi = {}
      if q_idx and q_idx.dof > 0 then
        for k = q_idx.start, q_idx.finish do
          table.insert(qi, q[k] or 0.0)
        end
      end

      local S_j = b.joint.subspace(qi) -- 6 x dof_j

      -- Transform subspace from joint j frame to body_id frame:
      -- We need transform X_{j -> body_id}
      -- X_world[j] transforms from 0 to j: v_j = X_0j * v_0
      -- X_world[body_id] transforms from 0 to body_id: v_body = X_0body * v_0
      -- Therefore: v_body = (X_0body * X_0j^(-1)) * v_j
      local X_0j_inv = matrix.inv(X_world[j])
      local X_j_to_body = X_world[body_id] * X_0j_inv

      local S_in_body = X_j_to_body * S_j

      if in_world_frame then
        -- Transform motion from body frame to world frame: X_0body^(-1)
        local X_0body_inv = matrix.inv(X_world[body_id])
        S_in_body = X_0body_inv * S_in_body
      end

      for col = 1, q_idx.dof do
        local col_idx = q_idx.start + col - 1
        for row = 1, 6 do
          J[row][col_idx] = S_in_body[row][col]
        end
      end
    end
  end

  return J
end

function jacobians.com_jacobian(m, q)
  -- 3 x dof Center of Mass Jacobian
  local _, X_world = kinematics.forward_kinematics(m, q)
  local J_com = matrix.zeros(3, m.dof)
  local total_mass = 0.0

  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    total_mass = total_mass + b.mass
    local J_body = jacobians.geometric_jacobian(m, q, i, true)

    -- Linear velocity rows in world frame: rows 4..6
    for col = 1, m.dof do
      for row = 1, 3 do
        J_com[row][col] = J_com[row][col] + b.mass * J_body[3 + row][col]
      end
    end
  end

  if total_mass > 0 then
    local inv_m = 1.0 / total_mass
    for row = 1, 3 do
      for col = 1, m.dof do
        J_com[row][col] = J_com[row][col] * inv_m
      end
    end
  end

  return J_com
end

function jacobians.nullspace_projector(J, tol)
  -- Computes Moore-Penrose pseudo-inverse J^+ and nullspace projector (I - J^+ * J) via SVD
  tol = tol or 1e-7
  local J_pinv = decomp.pinv(J, tol)
  local I = matrix.eye(J.cols)
  local N = I - (J_pinv * J)

  return J_pinv, N
end

function jacobians.inverse_kinematics_step(m, q, body_id, target_pos, target_R, step_size, damping)
  step_size = step_size or 0.2
  damping = damping or 1e-3

  local _, X_world = kinematics.forward_kinematics(m, q)
  local X_0b_inv = matrix.inv(X_world[body_id])
  local curr_pos = vec3.new(-X_0b_inv[5][3], X_0b_inv[4][3], -X_0b_inv[4][2])

  local delta_pos = target_pos - curr_pos
  local delta_x = { 0, 0, 0, delta_pos.x, delta_pos.y, delta_pos.z }

  local J = jacobians.geometric_jacobian(m, q, body_id, true)
  local J_pinv, _ = jacobians.nullspace_projector(J, damping)

  local delta_q = matrix.mul_vec(J_pinv, delta_x)
  local q_next = {}
  for i = 1, m.dof do
    q_next[i] = (q[i] or 0.0) + delta_q[i] * step_size
  end

  local err = vec3.norm(delta_pos)
  return q_next, err
end

return jacobians
