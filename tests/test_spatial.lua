local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")

runner.suite("Featherstone Spatial 6D Algebra")

runner.test("Spatial cross product duality relation", function()
  -- Identity: (v x m) . f = -m . (v x* f)
  local v = { 1.2, -0.5, 2.3, 4.1, -1.9, 0.8 }
  local m = { -0.7, 3.1, 1.4, -2.5, 0.4, 5.2 }
  local f = { 10.0, -15.0, 8.0, 20.0, -5.0, 12.0 }

  local v_cross_m = spatial.cross_motion(v, m)
  local lhs = spatial.dot(f, v_cross_m)

  local v_cross_star_f = spatial.cross_force(v, f)
  local rhs = -spatial.dot(v_cross_star_f, m)

  runner.assert_near(lhs, rhs, 1e-10)
end)

runner.test("Spatial coordinate transform composition", function()
  -- Frame A -> B with rotation R1 and offset p1
  -- Frame B -> C with rotation R2 and offset p2
  local R1 = mat3.from_axis_angle(vec3.unit_z(), 0.5)
  local p1 = vec3.new(1.0, 2.0, 0.0)
  local X_AB = spatial.transform(R1, p1)

  local R2 = mat3.from_axis_angle(vec3.unit_x(), 0.8)
  local p2 = vec3.new(0.0, -1.5, 3.0)
  local X_BC = spatial.transform(R2, p2)

  -- Compound transform X_AC = X_BC * X_AB (or X_AB * X_BC depending on convention)
  -- If v_A = X_AB * v_B and v_B = X_BC * v_C, then v_A = (X_AB * X_BC) * v_C
  local X_AC = X_AB * X_BC

  local v_C = { 0.2, -0.4, 1.1, 2.0, -3.0, 5.0 }
  local v_B = spatial.transform_motion(X_BC, v_C)
  local v_A1 = spatial.transform_motion(X_AB, v_B)
  local v_A2 = spatial.transform_motion(X_AC, v_C)

  runner.assert_vector_near(v_A1, v_A2, 1e-10)
end)

runner.test("Spatial inertia matrix and kinetic energy", function()
  local mass = 5.0
  local com = vec3.new(0.1, -0.2, 0.05)
  local I_rot = mat3.from_diag(0.5, 0.8, 0.4)

  local I_sp = spatial.inertia(mass, com, I_rot)
  runner.assert_equal(I_sp.rows, 6)
  runner.assert_equal(I_sp.cols, 6)

  -- Check symmetry: I_sp^T == I_sp
  for i = 1, 6 do
    for j = 1, 6 do
      runner.assert_near(I_sp[i][j], I_sp[j][i], 1e-12)
    end
  end

  -- Check positive kinetic energy for non-zero velocity
  local v = { 1.0, 2.0, -1.0, 0.5, -0.5, 3.0 }
  local T = spatial.kinetic_energy(I_sp, v)
  runner.assert_true(T > 0.0, "Kinetic energy must be strictly positive")
end)
