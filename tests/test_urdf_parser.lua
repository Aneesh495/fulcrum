local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")
local urdf_parser = require("fulcrum.models.urdf_parser")
local kinematics = require("fulcrum.dynamics.kinematics")
local rnea = require("fulcrum.dynamics.rnea")
local aba = require("fulcrum.dynamics.aba")

runner.suite("Robot Builder and URDF Parser")

runner.test("Fluent builder constructs 3-DOF articulated arm", function()
  local b = builder.new("arm_3dof")
  b:add_base("base_link", { mass = 0.0 })

  b:add_link({
    name = "link1",
    parent = "base_link",
    joint_type = "revolute",
    axis = vec3.unit_z(),
    origin = { xyz = vec3.new(0, 0, 0.2), rpy = vec3.zero() },
    mass = 2.0,
    com = vec3.new(0, 0, 0.1),
    inertia = { ixx = 0.05, iyy = 0.05, izz = 0.02 }
  })

  b:add_link({
    name = "link2",
    parent = "link1",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0, 0, 0.3), rpy = vec3.zero() },
    mass = 1.5,
    com = vec3.new(0, 0, 0.15),
    inertia = { ixx = 0.04, iyy = 0.04, izz = 0.01 }
  })

  b:add_link({
    name = "link3",
    parent = "link2",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0, 0, 0.3), rpy = vec3.zero() },
    mass = 1.0,
    com = vec3.new(0, 0, 0.1),
    inertia = { ixx = 0.02, iyy = 0.02, izz = 0.005 }
  })

  local model, meta = b:build()

  runner.assert_equal(model.num_bodies, 3, "Model has 3 active moving bodies")
  runner.assert_equal(model.dof, 3, "Model has 3 degrees of freedom")
  runner.assert_equal(meta.name, "arm_3dof", "Model metadata name matches")

  -- Forward kinematics check
  local q = { 0.0, 0.0, 0.0 }
  local X_tree, X_world = kinematics.forward_kinematics(model, q)
  runner.assert_equal(#X_world, 3, "Computed world transforms for all 3 links")

  -- Inverse dynamics check
  local qd = { 0.0, 0.0, 0.0 }
  local qdd = { 0.0, 0.0, 0.0 }
  local tau = rnea.inverse_dynamics(model, q, qd, qdd)
  runner.assert_equal(#tau, 3, "Torque vector has 3 elements")
end)

runner.test("URDF parser extracts links, joints, and constructs model", function()
  local urdf_sample = [[
<?xml version="1.0"?>
<robot name="scara_robot">
  <link name="base_link">
    <inertial>
      <mass value="0.0"/>
      <origin xyz="0 0 0" rpy="0 0 0"/>
    </inertial>
  </link>

  <link name="arm1">
    <inertial>
      <mass value="3.0"/>
      <origin xyz="0.2 0 0" rpy="0 0 0"/>
      <inertia ixx="0.01" ixy="0" ixz="0" iyy="0.1" iyz="0" izz="0.1"/>
    </inertial>
  </link>

  <link name="arm2">
    <inertial>
      <mass value="2.0"/>
      <origin xyz="0.15 0 0" rpy="0 0 0"/>
      <inertia ixx="0.008" ixy="0" ixz="0" iyy="0.05" iyz="0" izz="0.05"/>
    </inertial>
  </link>

  <joint name="joint1" type="revolute">
    <parent link="base_link"/>
    <child link="arm1"/>
    <origin xyz="0 0 0.1" rpy="0 0 0"/>
    <axis xyz="0 0 1"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="10"/>
  </joint>

  <joint name="joint2" type="revolute">
    <parent link="arm1"/>
    <child link="arm2"/>
    <origin xyz="0.4 0 0" rpy="0 0 0"/>
    <axis xyz="0 0 1"/>
    <limit lower="-2.5" upper="2.5" effort="80" velocity="10"/>
  </joint>
</robot>
]]

  local model, meta = urdf_parser.parse_string(urdf_sample)

  runner.assert_equal(meta.name, "scara_robot", "Robot name parsed")
  runner.assert_equal(model.num_bodies, 2, "SCARA model has 2 moving bodies")
  runner.assert_equal(model.dof, 2, "SCARA model has 2 DOF")

  -- Forward dynamics roundtrip
  local q = { 0.5, -0.5 }
  local qd = { 0.1, -0.2 }
  local tau = { 1.5, -0.8 }

  local qdd = aba.forward_dynamics(model, q, qd, tau)
  local tau_check = rnea.inverse_dynamics(model, q, qd, qdd)

  runner.assert_near(tau[1], tau_check[1], 1e-6, "ABA and RNEA consistency on parsed URDF joint 1")
  runner.assert_near(tau[2], tau_check[2], 1e-6, "ABA and RNEA consistency on parsed URDF joint 2")
end)
