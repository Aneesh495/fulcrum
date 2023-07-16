local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")

local aba = {}

function aba.forward_dynamics(m, q, q_dot, tau, f_ext)
  -- Featherstone Articulated Body Algorithm: computes q_ddot in O(N)
  local v = {}
  local c = {}
  local IA = {}
  local pA = {}
  local X_tree = {}
  local S_list = {}
  local U = {}
  local D_inv = {}
  local u = {}
  local a = {}
  local q_ddot = {}

  for i = 1, m.dof do
    q_ddot[i] = 0.0
  end

  local a_base = {
    m.gravity[1],
    m.gravity[2],
    m.gravity[3],
    -m.gravity[4],
    -m.gravity[5],
    -m.gravity[6]
  }

  -- Pass 1: Outward pass to compute kinematics and initialize articulated bodies
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

    local X_J = b.joint.transform(qi)
    local X_pi = X_J * b.X_T
    X_tree[i] = X_pi

    local S = b.joint.subspace(qi)
    S_list[i] = S

    local v_joint = matrix.mul_vec(S, qd_i)

    if b.parent_id == 0 then
      v[i] = v_joint
    else
      local v_from_parent = spatial.transform_motion(X_pi, v[b.parent_id])
      v[i] = spatial.add(v_from_parent, v_joint)
    end

    -- Joint bias acceleration c_i = v_i x (S_i * qd_i) + c_J
    c[i] = spatial.cross_motion(v[i], v_joint)

    -- Initialize IA and pA
    IA[i] = matrix.clone(b.I_spatial)
    local I_v = matrix.mul_vec(b.I_spatial, v[i])
    pA[i] = spatial.cross_force(v[i], I_v)

    if f_ext and f_ext[i] then
      pA[i] = spatial.sub(pA[i], f_ext[i])
    end
  end

  -- Pass 2: Inward pass to compute articulated body inertias and bias forces
  for i = m.num_bodies, 1, -1 do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local S = S_list[i]
    local dof_i = q_idx and q_idx.dof or 0

    if dof_i > 0 then
      -- U_i = IA_i * S_i (6 x dof_i)
      U[i] = IA[i] * S

      -- D_i = S_i^T * U_i (dof_i x dof_i)
      local S_t = matrix.transpose(S)
      local D = S_t * U[i]
      D_inv[i] = matrix.inv(D)

      -- u_i = tau_i - S_i^T * pA_i
      local St_pA = matrix.mul_vec(S_t, pA[i])
      u[i] = {}
      for k = 1, dof_i do
        local global_idx = q_idx.start + k - 1
        local tau_val = tau and tau[global_idx] or 0.0
        u[i][k] = tau_val - St_pA[k]
      end

      -- If body has parent, accumulate to parent
      if b.parent_id > 0 then
        -- Ia = IA_i - U_i * D_inv_i * U_i^T
        local Ut = matrix.transpose(U[i])
        local Ia = IA[i] - (U[i] * D_inv[i] * Ut)

        -- pa = pA_i + Ia * c_i + U_i * D_inv_i * u_i
        local Ia_c = matrix.mul_vec(Ia, c[i])
        local Dinv_u = matrix.mul_vec(D_inv[i], u[i])
        local U_Dinv_u = matrix.mul_vec(U[i], Dinv_u)
        local pa = spatial.add(pA[i], spatial.add(Ia_c, U_Dinv_u))

        -- Transform to parent frame: X_pi^T * Ia * X_pi, and X_pi^T * pa
        local X_pi_t = matrix.transpose(X_tree[i])
        local transformed_Ia = X_pi_t * Ia * X_tree[i]
        local transformed_pa = matrix.mul_vec(X_pi_t, pa)

        IA[b.parent_id] = IA[b.parent_id] + transformed_Ia
        pA[b.parent_id] = spatial.add(pA[b.parent_id], transformed_pa)
      end
    else
      -- Fixed joint
      if b.parent_id > 0 then
        local X_pi_t = matrix.transpose(X_tree[i])
        local transformed_Ia = X_pi_t * IA[i] * X_tree[i]
        local transformed_pa = matrix.mul_vec(X_pi_t, pA[i])
        IA[b.parent_id] = IA[b.parent_id] + transformed_Ia
        pA[b.parent_id] = spatial.add(pA[b.parent_id], transformed_pa)
      end
    end
  end

  -- Pass 3: Outward pass to compute joint accelerations
  for i = 1, m.num_bodies do
    local b = m.bodies[i]
    local q_idx = m.q_indices[i]
    local S = S_list[i]
    local dof_i = q_idx and q_idx.dof or 0

    local a_parent = (b.parent_id == 0) and a_base or a[b.parent_id]
    local a_prime = spatial.add(spatial.transform_motion(X_tree[i], a_parent), c[i])

    if dof_i > 0 then
      -- qdd_i = D_inv_i * (u_i - U_i^T * a_prime)
      local Ut = matrix.transpose(U[i])
      local Ut_aprime = matrix.mul_vec(Ut, a_prime)

      local rhs = {}
      for k = 1, dof_i do
        rhs[k] = u[i][k] - Ut_aprime[k]
      end

      local qdd_i = matrix.mul_vec(D_inv[i], rhs)
      for k = 1, dof_i do
        local global_idx = q_idx.start + k - 1
        q_ddot[global_idx] = qdd_i[k]
      end

      -- a_i = a_prime + S_i * qdd_i
      local S_qdd = matrix.mul_vec(S, qdd_i)
      a[i] = spatial.add(a_prime, S_qdd)
    else
      a[i] = a_prime
    end
  end

  return q_ddot
end

return aba
