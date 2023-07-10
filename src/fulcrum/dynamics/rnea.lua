local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")

local rnea = {}

function rnea.inverse_dynamics(m, q, q_dot, q_ddot, f_ext)
  -- Recursive Newton-Euler Algorithm: computes tau in O(N)
  -- tau = M(q) * q_ddot + C(q, q_dot)*q_dot + g(q) - J^T * f_ext
  local v = {}
  local a = {}
  local f = {}
  local X_tree = {}
  local S_list = {}
  local tau = {}

  -- Initialize tau vector
  for i = 1, m.dof do
    tau[i] = 0.0
  end

  -- Base acceleration initialized to -gravity (D'Alembert's principle)
  local a_base = {
    m.gravity[1],
    m.gravity[2],
    m.gravity[3],
    -m.gravity[4],
    -m.gravity[5],
    -m.gravity[6]
  }

  -- Outward pass (1 to N)
  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local qi = {}
    local qd_i = {}
    local qdd_i = {}

    if q_idx and q_idx.dof > 0 then
      for k = q_idx.start, q_idx.finish do
        table.insert(qi, q[k] or 0.0)
        table.insert(qd_i, q_dot[k] or 0.0)
        table.insert(qdd_i, q_ddot[k] or 0.0)
      end
    end

    local X_J = b.joint.transform(qi)
    local X_pi = X_J * b.X_T
    X_tree[i] = X_pi

    local S = b.joint.subspace(qi)
    S_list[i] = S

    local v_joint = matrix.mul_vec(S, qd_i)
    local a_joint = matrix.mul_vec(S, qdd_i)

    if b.parent_id == 0 then
      v[i] = v_joint
      local a_from_parent = spatial.transform_motion(X_pi, a_base)
      local v_cross_vJ = spatial.cross_motion(v[i], v_joint)
      a[i] = spatial.add(a_from_parent, spatial.add(a_joint, v_cross_vJ))
    else
      local v_from_parent = spatial.transform_motion(X_pi, v[b.parent_id])
      v[i] = spatial.add(v_from_parent, v_joint)

      local a_from_parent = spatial.transform_motion(X_pi, a[b.parent_id])
      local v_cross_vJ = spatial.cross_motion(v[i], v_joint)
      a[i] = spatial.add(a_from_parent, spatial.add(a_joint, v_cross_vJ))
    end

    -- Spatial force required to accelerate body: f_i = I_i * a_i + v_i x* (I_i * v_i)
    local I_a = matrix.mul_vec(b.I_spatial, a[i])
    local I_v = matrix.mul_vec(b.I_spatial, v[i])
    local v_cross_Iv = spatial.cross_force(v[i], I_v)
    f[i] = spatial.add(I_a, v_cross_Iv)

    -- Subtract external forces if specified on this body
    if f_ext and f_ext[i] then
      f[i] = spatial.sub(f[i], f_ext[i])
    end
  end

  -- Inward pass (N down to 1)
  for i = m.num_bodies, 1, -1 do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local S = S_list[i]

    -- Joint torque: tau_i = S_i^T * f_i
    if q_idx and q_idx.dof > 0 then
      for col = 1, q_idx.dof do
        local sum = 0.0
        for row = 1, 6 do
          sum = sum + S[row][col] * f[i][row]
        end
        local global_idx = q_idx.start + col - 1
        tau[global_idx] = sum
      end
    end

    -- Propagate force to parent: f_p = f_p + X_pi^* * f_i
    if b.parent_id > 0 then
      local f_to_parent = spatial.transform_force(X_tree[i], f[i])
      f[b.parent_id] = spatial.add(f[b.parent_id], f_to_parent)
    end
  end

  return tau
end

function rnea.gravity_terms(m, q)
  local zero = {}
  for i = 1, m.dof do zero[i] = 0.0 end
  return rnea.inverse_dynamics(m, q, zero, zero)
end

function rnea.coriolis_centrifugal(m, q, q_dot)
  local zero = {}
  for i = 1, m.dof do zero[i] = 0.0 end
  -- C(q, q_dot)*q_dot = RNEA(q, q_dot, 0) - g(q)
  local g = rnea.gravity_terms(m, q)
  local rnea_qd = rnea.inverse_dynamics(m, q, q_dot, zero)
  local c = {}
  for i = 1, m.dof do
    c[i] = rnea_qd[i] - g[i]
  end
  return c
end

return rnea
