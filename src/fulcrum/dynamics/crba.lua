local matrix = require("fulcrum.math.matrix")

local crba = {}

function crba.inertia_matrix(m, q)
  -- Composite Rigid Body Algorithm (CRBA) for M(q)
  local X_tree = {}
  local S_list = {}
  local Ic = {}

  -- Pass 1: Forward pass to compute transforms and initialize composite inertias
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
    local X_pi = X_J * b.X_T
    X_tree[i] = X_pi

    local S = b.joint.subspace(qi)
    S_list[i] = S

    -- Clone spatial inertia
    Ic[i] = matrix.clone(b.I_spatial)
  end

  -- Pass 2: Inward pass to accumulate composite subtree inertias
  for i = m.num_bodies, 1, -1 do
    local b = m.bodies[i]
    if b.parent_id > 0 then
      -- Ic[parent] = Ic[parent] + X_pi^T * Ic[i] * X_pi
      local X_pi_t = matrix.transpose(X_tree[i])
      local transformed_Ic = X_pi_t * Ic[i] * X_tree[i]
      Ic[b.parent_id] = Ic[b.parent_id] + transformed_Ic
    end
  end

  -- Pass 3: Assemble joint-space inertia matrix M
  local M = matrix.zeros(m.dof, m.dof)

  for i = 1, m.num_bodies do
    local q_idx_i = m.q_indices[i]
    if q_idx_i and q_idx_i.dof > 0 then
      local S_i = S_list[i]
      -- F = Ic[i] * S_i (6 x dof_i)
      local F = Ic[i] * S_i

      -- Diagonal block M_ii = S_i^T * F
      local S_i_t = matrix.transpose(S_i)
      local M_ii = S_i_t * F
      for r = 1, q_idx_i.dof do
        for c = 1, q_idx_i.dof do
          local gi = q_idx_i.start + r - 1
          local gj = q_idx_i.start + c - 1
          M[gi][gj] = M_ii[r][c]
        end
      end

      -- Off-diagonal blocks: traverse up ancestor tree of body i
      local curr_j = m.bodies[i].parent_id
      local F_up = F
      local curr_child = i

      while curr_j > 0 do
        local q_idx_j = m.q_indices[curr_j]
        -- F_up = X_tree[curr_child]^T * F_up
        local X_child_t = matrix.transpose(X_tree[curr_child])
        F_up = X_child_t * F_up

        if q_idx_j and q_idx_j.dof > 0 then
          local S_j = S_list[curr_j]
          local S_j_t = matrix.transpose(S_j)
          local M_ji = S_j_t * F_up

          for r = 1, q_idx_j.dof do
            for c = 1, q_idx_i.dof do
              local gj = q_idx_j.start + r - 1
              local gi = q_idx_i.start + c - 1
              M[gj][gi] = M_ji[r][c]
              M[gi][gj] = M_ji[r][c] -- Symmetry
            end
          end
        end

        curr_child = curr_j
        curr_j = m.bodies[curr_j].parent_id
      end
    end
  end

  return M
end

crba.mass_matrix = crba.inertia_matrix

return crba
