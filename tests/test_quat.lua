local runner = require("tests.test_framework")
local quat = require("fulcrum.math.quat")
local vec3 = require("fulcrum.math.vec3")

runner.suite("Quaternion Algebra")

runner.test("Quaternion identity and norm", function()
  local q = quat.identity()
  runner.assert_equal(q.w, 1.0)
  runner.assert_equal(q.x, 0.0)
  runner.assert_equal(q.y, 0.0)
  runner.assert_equal(q.z, 0.0)
  runner.assert_near(quat.norm(q), 1.0)
end)

runner.test("Hamilton product multiplication", function()
  -- i * j = k
  local qi = quat.new(0, 1, 0, 0)
  local qj = quat.new(0, 0, 1, 0)
  local qk = quat.new(0, 0, 0, 1)

  local prod = qi * qj
  runner.assert_near(prod.w, qk.w)
  runner.assert_near(prod.x, qk.x)
  runner.assert_near(prod.y, qk.y)
  runner.assert_near(prod.z, qk.z)

  -- j * k = i
  local jk = qj * qk
  runner.assert_near(jk.w, qi.w)
  runner.assert_near(jk.x, qi.x)
  runner.assert_near(jk.y, qi.y)
  runner.assert_near(jk.z, qi.z)

  -- i^2 = j^2 = k^2 = ijk = -1
  local i2 = qi * qi
  runner.assert_near(i2.w, -1.0)
  runner.assert_near(i2.x, 0.0)
  runner.assert_near(i2.y, 0.0)
  runner.assert_near(i2.z, 0.0)
end)

runner.test("Axis-angle roundtrip", function()
  local axis = vec3.unit(vec3.new(1.0, 2.0, 3.0))
  local angle = 1.25 -- radians

  local q = quat.from_axis_angle(axis, angle)
  local rec_axis, rec_angle = quat.to_axis_angle(q)

  runner.assert_near(rec_angle, angle, 1e-10)
  runner.assert_near(rec_axis.x, axis.x, 1e-10)
  runner.assert_near(rec_axis.y, axis.y, 1e-10)
  runner.assert_near(rec_axis.z, axis.z, 1e-10)
end)

runner.test("Vector rotation", function()
  -- Rotate unit X by 90 degrees around Z, should yield unit Y
  local q = quat.from_axis_angle(vec3.unit_z(), math.pi * 0.5)
  local v = vec3.unit_x()
  local rotated = quat.rotate_vector(q, v)

  runner.assert_near(rotated.x, 0.0, 1e-10)
  runner.assert_near(rotated.y, 1.0, 1e-10)
  runner.assert_near(rotated.z, 0.0, 1e-10)
end)

runner.test("Rotation matrix roundtrip", function()
  local axis = vec3.unit(vec3.new(0.5, -0.8, 1.2))
  local angle = 0.85
  local q = quat.from_axis_angle(axis, angle)

  local R = quat.to_rotation_matrix(q)
  local q_recovered = quat.from_rotation_matrix(R)

  -- Quaternion equality accounts for antipodal sign
  local dot = math.abs(q.w * q_recovered.w + q.x * q_recovered.x + q.y * q_recovered.y + q.z * q_recovered.z)
  runner.assert_near(dot, 1.0, 1e-7)
end)

runner.test("Slerp interpolation", function()
  local q1 = quat.identity()
  local q2 = quat.from_axis_angle(vec3.unit_y(), math.pi * 0.5)

  local q_mid = quat.slerp(q1, q2, 0.5)
  local axis, angle = quat.to_axis_angle(q_mid)

  runner.assert_near(angle, math.pi * 0.25, 1e-7)
  runner.assert_near(axis.y, 1.0, 1e-7)
end)
