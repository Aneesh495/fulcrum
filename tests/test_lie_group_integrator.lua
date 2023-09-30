local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local so3 = require("fulcrum.math.so3")
local se3 = require("fulcrum.math.se3")
local lgvi = require("fulcrum.integrators.lie_group_integrator")

runner.suite("Lie Group Variational Integrator (LGVI)")

runner.test("Energy and angular momentum conservation on SO(3)", function()
  -- Asymmetric top (J1 < J2 < J3)
  local J_diag = vec3.new(1.0, 2.0, 3.0)
  local R = mat3.identity()
  -- Initial angular velocity near intermediate axis (demonstrates tennis racket effect)
  local Omega = vec3.new(0.01, 2.0, 0.01)

  local kinetic_energy = function(w)
    return 0.5 * (J_diag.x * w.x * w.x + J_diag.y * w.y * w.y + J_diag.z * w.z * w.z)
  end

  local spatial_momentum_norm = function(rot, w)
    local Pi_body = vec3.new(J_diag.x * w.x, J_diag.y * w.y, J_diag.z * w.z)
    local Pi_spatial = rot * Pi_body
    return vec3.norm(Pi_spatial)
  end

  local E0 = kinetic_energy(Omega)
  local L0 = spatial_momentum_norm(R, Omega)
  local dt = 0.005

  for step = 1, 500 do
    R, Omega = lgvi.step_so3(R, Omega, J_diag, nil, dt)
  end

  local E_final = kinetic_energy(Omega)
  local L_final = spatial_momentum_norm(R, Omega)

  -- Energy conserved within symplectic fluctuation bound
  runner.assert_near(E_final, E0, 1e-3)
  -- Spatial angular momentum magnitude is exactly conserved (Noether theorem)
  runner.assert_near(L_final, L0, 1e-8)

  -- Check R remains exactly in SO(3)
  local RtR = mat3.transpose(R) * R
  local I = mat3.identity()
  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(RtR[i][j], I[i][j], 1e-10)
    end
  end
  runner.assert_near(mat3.det(R), 1.0, 1e-10)
end)

runner.test("SE(3) rigid body free flight", function()
  local mass = 10.0
  local J_diag = vec3.new(2.0, 2.0, 2.0)
  local T = se3.identity()
  local twist = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 }
  local dt = 0.01

  for step = 1, 100 do
    T, twist = lgvi.step_se3(T, twist, mass, J_diag, nil, dt)
  end

  -- Body with vx=1.0 and wz=1.0 follows circular arc: x(1) = sin(1) ~ 0.8415, y(1) = 1 - cos(1) ~ 0.4597
  runner.assert_near(T.p.x, math.sin(1.0), 0.01)
  runner.assert_near(T.p.y, 1.0 - math.cos(1.0), 0.01)
end)
