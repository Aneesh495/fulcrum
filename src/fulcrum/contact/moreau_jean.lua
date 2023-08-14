local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")
local rnea = require("fulcrum.dynamics.rnea")
local crba = require("fulcrum.dynamics.crba")
local lemke = require("fulcrum.contact.lemke")
local pgs = require("fulcrum.contact.pgs")

local moreau_jean = {}

function moreau_jean.step(m, q, q_dot, tau, dt, contacts, restitution, mu, solver_type)
  -- Moreau-Jean velocity-impulse time-stepping scheme
  -- contacts: array of contact objects { body_id, normal, J_normal, J_tangent1, J_tangent2, gap }
  solver_type = solver_type or "lemke"
  restitution = restitution or 0.0
  mu = mu or 0.5
  local dof = m.dof

  local M = crba.inertia_matrix(m, q)
  local c = rnea.coriolis_centrifugal(m, q, q_dot)
  local grav = rnea.gravity_terms(m, q)

  -- Free velocity v_free (velocity without contact impulses):
  -- M * (v_free - v_k) = dt * (tau - c - grav)
  local tau_smooth = {}
  for i = 1, dof do
    local applied = tau and tau[i] or 0.0
    tau_smooth[i] = applied - c[i] - grav[i]
  end

  local M_inv = matrix.inv(M)
  local dv_free = matrix.mul_vec(M_inv, tau_smooth)
  local v_free = {}
  for i = 1, dof do
    v_free[i] = (q_dot[i] or 0.0) + dv_free[i] * dt
  end

  if not contacts or #contacts == 0 then
    -- No active contacts: pure smooth update
    local q_next = {}
    for i = 1, dof do
      q_next[i] = (q[i] or 0.0) + v_free[i] * dt
    end
    return q_next, v_free, {}
  end

  -- Assemble contact Jacobian J_c: rows for normal contacts
  local num_contacts = #contacts
  local J_n = matrix.zeros(num_contacts, dof)
  local q_lcp = {}

  for i = 1, num_contacts do
    local ct = contacts[i]
    for j = 1, dof do
      J_n[i][j] = ct.J_normal[j] or 0.0
    end

    -- Normal pre-impact velocity: vn_k = J_n[i] * q_dot
    local vn_k = 0.0
    local vn_free = 0.0
    for j = 1, dof do
      vn_k = vn_k + J_n[i][j] * (q_dot[j] or 0.0)
      vn_free = vn_free + J_n[i][j] * v_free[j]
    end

    -- Signorini velocity offset: q_lcp[i] = vn_free - e * |vn_k|
    q_lcp[i] = vn_free - restitution * math.max(0.0, -vn_k)
  end

  -- Delassus matrix: W = J_n * M^(-1) * J_n^T
  local J_n_t = matrix.transpose(J_n)
  local W = J_n * M_inv * J_n_t

  -- Regularize Delassus matrix slightly
  for i = 1, num_contacts do
    W[i][i] = W[i][i] + 1e-7
  end

  -- Solve LCP: w = W * P_n + q_lcp >= 0, P_n >= 0, w^T * P_n = 0
  local status, P_n, w, _ = lemke.solve(W, q_lcp)

  -- If Lemke encounters ray or max iter, fallback to PGS
  if status ~= lemke.STATUS_SOLVED or not P_n then
    P_n = pgs.solve_box(W, q_lcp, nil, nil, 100, 1e-6, 1.0)
  end

  -- Post-impact velocity: v_{k+1} = v_free + M^(-1) * J_n^T * P_n
  local Jn_t_Pn = matrix.mul_vec(J_n_t, P_n)
  local dv_contact = matrix.mul_vec(M_inv, Jn_t_Pn)

  local v_next = {}
  local q_next = {}
  for i = 1, dof do
    v_next[i] = v_free[i] + dv_contact[i]
    q_next[i] = (q[i] or 0.0) + v_next[i] * dt
  end

  return q_next, v_next, P_n
end

return moreau_jean
