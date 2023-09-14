local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")

local manipulator = {}

function manipulator.create(params)
  params = params or {}
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("manipulator_7dof", gravity)
  b:add_base("base_link", { mass = 0.0 })

  -- Link 1: Shoulder Yaw
  b:add_link({
    name = "link1",
    parent = "base_link",
    joint_name = "joint1_yaw",
    joint_type = "revolute",
    axis = vec3.unit_z(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.35), rpy = vec3.zero() },
    mass = 4.0,
    com = vec3.new(0.0, 0.0, 0.1),
    inertia = { ixx = 0.08, iyy = 0.08, izz = 0.04 }
  })

  -- Link 2: Shoulder Pitch
  b:add_link({
    name = "link2",
    parent = "link1",
    joint_name = "joint2_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.2), rpy = vec3.zero() },
    mass = 4.0,
    com = vec3.new(0.0, 0.0, 0.15),
    inertia = { ixx = 0.09, iyy = 0.09, izz = 0.03 }
  })

  -- Link 3: Shoulder Roll
  b:add_link({
    name = "link3",
    parent = "link2",
    joint_name = "joint3_roll",
    joint_type = "revolute",
    axis = vec3.unit_z(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.3), rpy = vec3.zero() },
    mass = 3.0,
    com = vec3.new(0.0, 0.0, 0.1),
    inertia = { ixx = 0.05, iyy = 0.05, izz = 0.02 }
  })

  -- Link 4: Elbow Pitch
  b:add_link({
    name = "link4",
    parent = "link3",
    joint_name = "joint4_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.2), rpy = vec3.zero() },
    mass = 2.5,
    com = vec3.new(0.0, 0.0, 0.1),
    inertia = { ixx = 0.04, iyy = 0.04, izz = 0.015 }
  })

  -- Link 5: Wrist Roll
  b:add_link({
    name = "link5",
    parent = "link4",
    joint_name = "joint5_roll",
    joint_type = "revolute",
    axis = vec3.unit_z(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.25), rpy = vec3.zero() },
    mass = 2.0,
    com = vec3.new(0.0, 0.0, 0.08),
    inertia = { ixx = 0.025, iyy = 0.025, izz = 0.01 }
  })

  -- Link 6: Wrist Pitch
  b:add_link({
    name = "link6",
    parent = "link5",
    joint_name = "joint6_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.15), rpy = vec3.zero() },
    mass = 1.5,
    com = vec3.new(0.0, 0.0, 0.05),
    inertia = { ixx = 0.015, iyy = 0.015, izz = 0.008 }
  })

  -- Link 7: Wrist Yaw / Tool Flange
  b:add_link({
    name = "link7",
    parent = "link6",
    joint_name = "joint7_yaw",
    joint_type = "revolute",
    axis = vec3.unit_z(),
    origin = { xyz = vec3.new(0.0, 0.0, 0.12), rpy = vec3.zero() },
    mass = 0.8,
    com = vec3.new(0.0, 0.0, 0.03),
    inertia = { ixx = 0.006, iyy = 0.006, izz = 0.004 }
  })

  local model, meta = b:build()
  meta.dof = 7
  meta.end_effector_id = model.num_bodies
  return model, meta
end

return manipulator
