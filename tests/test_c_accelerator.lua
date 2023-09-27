local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local spatial = require("fulcrum.math.spatial")
local c_accel = require("fulcrum.ffi.c_accelerator")

runner.suite("High-Performance C Accelerator & FFI Kernels")

runner.test("Vector3 and Matrix3 C kernels consistency", function()
  local a = vec3.new(1.2, -3.4, 5.6)
  local b = vec3.new(-0.7, 2.1, -4.3)

  local cross_accel = c_accel.vec3_cross(a, b)
  local cross_pure = vec3.cross(a, b)

  runner.assert_near(cross_accel.x, cross_pure.x, 1e-9, "Cross product X matches")
  runner.assert_near(cross_accel.y, cross_pure.y, 1e-9, "Cross product Y matches")
  runner.assert_near(cross_accel.z, cross_pure.z, 1e-9, "Cross product Z matches")

  local dot_accel = c_accel.vec3_dot(a, b)
  local dot_pure = vec3.dot(a, b)
  runner.assert_near(dot_accel, dot_pure, 1e-9, "Dot product matches")

  local A = mat3.from_euler_xyz(0.2, -0.4, 0.6)
  local B = mat3.from_euler_xyz(-0.1, 0.3, -0.5)
  local mul_accel = c_accel.mat3_mul(A, B)
  local mul_pure = mat3.mul(A, B)

  for i = 1, 3 do
    for j = 1, 3 do
      runner.assert_near(mul_accel[i][j], mul_pure[i][j], 1e-9, "Matrix mul elements match")
    end
  end
end)

runner.test("Spatial 6D Featherstone transforms and kinetic energy", function()
  local R = mat3.from_axis_angle(vec3.unit_z(), 0.5)
  local r = vec3.new(1.0, 2.0, 3.0)
  local X = spatial.transform(R, r)
  local v = { 0.1, -0.2, 0.3, 1.0, -1.5, 2.0 }

  local v_accel = c_accel.spatial_transform(X, v)
  local v_pure = spatial.transform_motion(X, v)

  for i = 1, 6 do
    runner.assert_near(v_accel[i], v_pure[i], 1e-9, "Motion transform component " .. i .. " matches")
  end

  local f = { 10.0, -20.0, 30.0, 5.0, -2.0, 1.0 }
  local f_accel = c_accel.spatial_transform_force(X, f)
  local f_pure = spatial.transform_force(X, f)

  for i = 1, 6 do
    runner.assert_near(f_accel[i], f_pure[i], 1e-9, "Force transform component " .. i .. " matches")
  end

  local I_sp = spatial.inertia(2.5, vec3.new(0.1, 0.2, 0.3), mat3.from_diag(0.5, 0.6, 0.7))
  local T_accel = c_accel.spatial_kinetic_energy(I_sp, v)
  local T_pure = spatial.kinetic_energy(I_sp, v)
  runner.assert_near(T_accel, T_pure, 1e-9, "Spatial kinetic energy matches")
end)

runner.test("Projected Gauss-Seidel LCP contact solve", function()
  -- Delassus matrix for contact system
  local A = {
    { 2.0, -1.0 },
    { -1.0, 2.0 }
  }
  local b = { -3.0, -1.0 }
  local lo = { 0.0, 0.0 }
  local hi = { 10.0, 10.0 }

  local sol, iters = c_accel.pgs_solve(A, b, lo, hi, { 0.0, 0.0 }, 100, 1e-6)
  runner.assert_true(iters < 50, "PGS converged in fewer than 50 iterations")

  -- Unconstrained minimum would be A x + b = 0 => x = A^-1 (-b)
  -- A^-1 = 1/3 * [2, 1; 1, 2], -b = [3, 1]^T
  -- x1 = 1/3 * (6 + 1) = 7/3 = 2.33333
  -- x2 = 1/3 * (3 + 2) = 5/3 = 1.66667
  runner.assert_near(sol[1], 7.0 / 3.0, 1e-4, "PGS solved x1")
  runner.assert_near(sol[2], 5.0 / 3.0, 1e-4, "PGS solved x2")
end)
