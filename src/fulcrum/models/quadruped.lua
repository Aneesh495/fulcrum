local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")

local quadruped = {}

function quadruped.create(params)
  params = params or {}
  local is_floating = (params.floating ~= false)
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local trunk_len = params.trunk_len or 0.4
  local trunk_width = params.trunk_width or 0.2
  local trunk_height = params.trunk_height or 0.1
  local trunk_mass = params.trunk_mass or 6.0

  local thigh_len = params.thigh_len or 0.2
  local shin_len = params.shin_len or 0.2

  local b = builder.new("quadruped", gravity)
  b:add_base("world", { mass = 0.0 })

  -- Trunk
  local I_xx = (1.0 / 12.0) * trunk_mass * (trunk_width * trunk_width + trunk_height * trunk_height)
  local I_yy = (1.0 / 12.0) * trunk_mass * (trunk_len * trunk_len + trunk_height * trunk_height)
  local I_zz = (1.0 / 12.0) * trunk_mass * (trunk_len * trunk_len + trunk_width * trunk_width)

  b:add_link({
    name = "trunk",
    parent = "world",
    joint_name = "floating_base_joint",
    joint_type = is_floating and "floating_base" or "fixed",
    origin = { xyz = vec3.new(0.0, 0.0, 0.4), rpy = vec3.zero() },
    mass = trunk_mass,
    com = vec3.zero(),
    inertia = { ixx = I_xx, iyy = I_yy, izz = I_zz }
  })

  -- Leg configurations: { name, x_offset, y_offset }
  local half_l = trunk_len * 0.5
  local half_w = trunk_width * 0.5
  local legs = {
    { prefix = "fl", x = half_l,  y = half_w },
    { prefix = "fr", x = half_l,  y = -half_w },
    { prefix = "rl", x = -half_l, y = half_w },
    { prefix = "rr", x = -half_l, y = -half_w }
  }

  for _, leg in ipairs(legs) do
    local p = leg.prefix
    -- Hip Roll
    b:add_link({
      name = p .. "_hip",
      parent = "trunk",
      joint_name = p .. "_hip_roll",
      joint_type = "revolute",
      axis = vec3.unit_x(),
      origin = { xyz = vec3.new(leg.x, leg.y, 0.0), rpy = vec3.zero() },
      mass = 0.8,
      com = vec3.zero(),
      inertia = { ixx = 0.002, iyy = 0.002, izz = 0.002 }
    })

    -- Thigh Pitch
    b:add_link({
      name = p .. "_thigh",
      parent = p .. "_hip",
      joint_name = p .. "_thigh_pitch",
      joint_type = "revolute",
      axis = vec3.unit_y(),
      origin = { xyz = vec3.new(0.0, 0.04 * (leg.y > 0 and 1 or -1), 0.0), rpy = vec3.zero() },
      mass = 1.0,
      com = vec3.new(0.0, 0.0, -0.5 * thigh_len),
      inertia = { ixx = 0.004, iyy = 0.004, izz = 0.001 }
    })

    -- Calf / Knee Pitch
    b:add_link({
      name = p .. "_calf",
      parent = p .. "_thigh",
      joint_name = p .. "_knee_pitch",
      joint_type = "revolute",
      axis = vec3.unit_y(),
      origin = { xyz = vec3.new(0.0, 0.0, -thigh_len), rpy = vec3.zero() },
      mass = 0.4,
      com = vec3.new(0.0, 0.0, -0.5 * shin_len),
      inertia = { ixx = 0.002, iyy = 0.002, izz = 0.0005 }
    })
  end

  local model, meta = b:build()
  meta.foot_ids = {
    fl = meta.link_to_id["fl_calf"],
    fr = meta.link_to_id["fr_calf"],
    rl = meta.link_to_id["rl_calf"],
    rr = meta.link_to_id["rr_calf"]
  }
  return model, meta
end

return quadruped
