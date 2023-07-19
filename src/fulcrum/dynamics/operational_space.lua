local matrix = require("fulcrum.math.matrix")
local decomp = require("fulcrum.math.decomp")
local rnea = require("fulcrum.dynamics.rnea")
local crba = require("fulcrum.dynamics.crba")

local op_space = {}

function op_space.inertia(M, J)
  -- Lambda(x) = (J * M^(-1) * J^T)^(-1)
  local M_inv = matrix.inv(M)
  local Jt = matrix.transpose(J)
  local J_Minv_Jt = J * M_inv * Jt

  -- Regularized pseudo-inverse or inverse of operational space mobility
  local Lambda = decomp.pinv(J_Minv_Jt, 1e-8)
  return Lambda, M_inv
end

function op_space.dynamically_consistent_inverse(M_inv, J, Lambda)
  -- J^# = M^(-1) * J^T * Lambda
  local Jt = matrix.transpose(J)
  local J_sharp = M_inv * Jt * Lambda
  return J_sharp
end

function op_space.nullspace_projector(J, J_sharp)
  -- N_T = I - J^T * (J^#)^T
  local n = J.cols
  local I = matrix.eye(n)
  local Jt = matrix.transpose(J)
  local J_sharp_t = matrix.transpose(J_sharp)
  local N_T = I - (Jt * J_sharp_t)
  return N_T
end

function op_space.task_torques(m, q, q_dot, J, F_task, tau_null)
  -- tau = J^T * F_task + (I - J^T * J^#T) * tau_null + g(q)
  local M = crba.inertia_matrix(m, q)
  local Lambda, M_inv = op_space.inertia(M, J)
  local J_sharp = op_space.dynamically_consistent_inverse(M_inv, J, Lambda)
  local N_T = op_space.nullspace_projector(J, J_sharp)

  local Jt = matrix.transpose(J)
  local tau_task = matrix.mul_vec(Jt, F_task)

  local tau_posture = { 0, 0, 0 }
  if tau_null then
    tau_posture = matrix.mul_vec(N_T, tau_null)
  else
    for i = 1, m.dof do tau_posture[i] = 0.0 end
  end

  local g = rnea.gravity_terms(m, q)

  local tau = {}
  for i = 1, m.dof do
    tau[i] = tau_task[i] + tau_posture[i] + g[i]
  end

  return tau, Lambda, J_sharp
end

return op_space
