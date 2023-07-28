local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local so3 = require("fulcrum.math.so3")
local se3 = require("fulcrum.math.se3")

local lgvi = {}

function lgvi.step_so3(R_k, Omega_k, J_diag, tau_k, dt, max_iter, tol)
  -- Lie Group Variational Integrator on SO(3) (Lee, Leok, McClamroch)
  -- R_k: 3x3 rotation matrix
  -- Omega_k: body angular velocity (vec3)
  -- J_diag: body principal moments of inertia [J1, J2, J3] (vec3)
  -- tau_k: external torque in body frame (vec3)
  max_iter = max_iter or 15
  tol = tol or 1e-12

  local J = mat3.from_diag(J_diag.x, J_diag.y, J_diag.z)
  local Pi_k = vec3.new(
    J_diag.x * Omega_k.x,
    J_diag.y * Omega_k.y,
    J_diag.z * Omega_k.z
  )

  local tau = tau_k or vec3.zero()
  local M_vec = Pi_k + (tau * dt)

  -- Initial guess for f_k = dt * Omega_k
  local f = Omega_k * dt

  -- Solve implicit discrete Euler-Poincaré equation:
  -- res(f) = (F * J - J * F^T)^vee - h * M_vec = 0
  for _ = 1, max_iter do
    local F = so3.exp(f)
    local FJ = F * J
    local JFt = J * mat3.transpose(F)
    local skew_diff = FJ - JFt

    local v_diff = so3.vee(skew_diff)
    local res = v_diff - (M_vec * dt)

    if vec3.norm(res) < tol then break end

    -- Analytical Jacobian of residual with respect to f:
    -- D_res = tr(F*J)*I - (F*J) (approximate or exact)
    local tr_FJ = mat3.trace(FJ)
    local Grad = mat3.new(
      tr_FJ - FJ[1][1], -FJ[1][2],        -FJ[1][3],
      -FJ[2][1],        tr_FJ - FJ[2][2], -FJ[2][3],
      -FJ[3][1],        -FJ[3][2],        tr_FJ - FJ[3][3]
    )

    local delta_f = mat3.inv(Grad) * res
    f = f - delta_f
  end

  local F_k = so3.exp(f)
  local R_next = R_k * F_k
  -- Re-orthonormalize strictly to prevent machine precision drift
  R_next = mat3.orthonormalize(R_next)

  -- Update angular momentum: Pi_{k+1} = F_k^T * (Pi_k + dt * tau)
  local Ft = mat3.transpose(F_k)
  local Pi_next = Ft * M_vec

  local Omega_next = vec3.new(
    Pi_next.x / J_diag.x,
    Pi_next.y / J_diag.y,
    Pi_next.z / J_diag.z
  )

  return R_next, Omega_next
end

function lgvi.step_se3(T_k, twist_k, mass, J_diag, wrench_k, dt)
  -- Full SE(3) Lie Group Variational Integrator
  -- T_k: SE(3) pose
  -- twist_k: [vx, vy, vz, wx, wy, wz] (linear and angular body velocities)
  -- mass: scalar mass
  -- J_diag: vec3 rotational inertias
  -- wrench_k: [Fx, Fy, Fz, Tx, Ty, Tz]
  wrench_k = wrench_k or { 0, 0, 0, 0, 0, 0 }

  local v_k = vec3.new(twist_k[1], twist_k[2], twist_k[3])
  local w_k = vec3.new(twist_k[4], twist_k[5], twist_k[6])

  local F_ext = vec3.new(wrench_k[1], wrench_k[2], wrench_k[3])
  local T_ext = vec3.new(wrench_k[4], wrench_k[5], wrench_k[6])

  -- Step rotation on SO(3)
  local R_next, w_next = lgvi.step_so3(T_k.R, w_k, J_diag, T_ext, dt)

  -- Step translation (linear momentum update)
  local p_lin = v_k * mass + (F_ext * dt)
  local v_next = p_lin / mass

  -- Position step using exponential screw
  local p_next = T_k.p + (T_k.R * (v_k * dt))

  local T_next = se3.new(R_next, p_next)
  local twist_next = {
    v_next.x, v_next.y, v_next.z,
    w_next.x, w_next.y, w_next.z
  }

  return T_next, twist_next
end

return lgvi
