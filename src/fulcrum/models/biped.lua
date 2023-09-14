local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")

local biped = {}

function biped.create(params)
  params = params or {}
  local is_floating = (params.floating ~= false)
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("biped_walker", gravity)
  b:add_base("world", { mass = 0.0 })

  -- Pelvis / Torso
  b:add_link({
    name = "torso",
    parent = "world",
    joint_name = "base_joint",
    joint_type = is_floating and "floating_base" or "fixed",
    origin = { xyz = vec3.new(0.0, 0.0, 0.8), rpy = vec3.zero() },
    mass = 12.0,
    com = vec3.new(0.0, 0.0, 0.15),
    inertia = { ixx = 0.3, iyy = 0.3, izz = 0.2 }
  })

  -- Left Leg: Hip, Knee, Ankle
  local leg_width = 0.12
  local thigh_len = 0.35
  local shin_len = 0.35
  local foot_len = 0.15

  b:add_link({
    name = "left_thigh",
    parent = "torso",
    joint_name = "left_hip_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, leg_width, 0.0), rpy = vec3.zero() },
    mass = 3.0,
    com = vec3.new(0.0, 0.0, -0.5 * thigh_len),
    inertia = { ixx = 0.04, iyy = 0.04, izz = 0.01 }
  })

  b:add_link({
    name = "left_shin",
    parent = "left_thigh",
    joint_name = "left_knee_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, -thigh_len), rpy = vec3.zero() },
    mass = 2.0,
    com = vec3.new(0.0, 0.0, -0.5 * shin_len),
    inertia = { ixx = 0.03, iyy = 0.03, izz = 0.008 }
  })

  b:add_link({
    name = "left_foot",
    parent = "left_shin",
    joint_name = "left_ankle_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, -shin_len), rpy = vec3.zero() },
    mass = 0.8,
    com = vec3.new(0.04, 0.0, -0.03),
    inertia = { ixx = 0.005, iyy = 0.005, izz = 0.005 }
  })

  -- Right Leg: Hip, Knee, Ankle
  b:add_link({
    name = "right_thigh",
    parent = "torso",
    joint_name = "right_hip_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, -leg_width, 0.0), rpy = vec3.zero() },
    mass = 3.0,
    com = vec3.new(0.0, 0.0, -0.5 * thigh_len),
    inertia = { ixx = 0.04, iyy = 0.04, izz = 0.01 }
  })

  b:add_link({
    name = "right_shin",
    parent = "right_thigh",
    joint_name = "right_knee_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, -thigh_len), rpy = vec3.zero() },
    mass = 2.0,
    com = vec3.new(0.0, 0.0, -0.5 * shin_len),
    inertia = { ixx = 0.03, iyy = 0.03, izz = 0.008 }
  })

  b:add_link({
    name = "right_foot",
    parent = "right_shin",
    joint_name = "right_ankle_pitch",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, -shin_len), rpy = vec3.zero() },
    mass = 0.8,
    com = vec3.new(0.04, 0.0, -0.03),
    inertia = { ixx = 0.005, iyy = 0.005, izz = 0.005 }
  })

  local model, meta = b:build()
  meta.left_foot_id = meta.link_to_id["left_foot"]
  meta.right_foot_id = meta.link_to_id["right_foot"]
  return model, meta
end

return biped
