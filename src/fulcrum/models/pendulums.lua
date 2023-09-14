local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")

local pendulums = {}

function pendulums.double_pendulum(params)
  params = params or {}
  local m1 = params.m1 or 1.0
  local m2 = params.m2 or 1.0
  local l1 = params.l1 or 1.0
  local l2 = params.l2 or 1.0
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("double_pendulum", gravity)
  b:add_base("world", { mass = 0.0 })

  local I1 = (1.0 / 12.0) * m1 * l1 * l1
  b:add_link({
    name = "link1",
    parent = "world",
    joint_name = "joint1",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.zero(), rpy = vec3.zero() },
    mass = m1,
    com = vec3.new(0.0, 0.0, -0.5 * l1),
    inertia = { ixx = I1, iyy = I1, izz = 0.001 }
  })

  local I2 = (1.0 / 12.0) * m2 * l2 * l2
  b:add_link({
    name = "link2",
    parent = "link1",
    joint_name = "joint2",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.new(0.0, 0.0, -l1), rpy = vec3.zero() },
    mass = m2,
    com = vec3.new(0.0, 0.0, -0.5 * l2),
    inertia = { ixx = I2, iyy = I2, izz = 0.001 }
  })

  local model, meta = b:build()
  meta.m1, meta.m2 = m1, m2
  meta.l1, meta.l2 = l1, l2
  return model, meta
end

function pendulums.triple_pendulum(params)
  params = params or {}
  local m = params.masses or { 1.0, 1.0, 1.0 }
  local l = params.lengths or { 1.0, 1.0, 1.0 }
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("triple_pendulum", gravity)
  b:add_base("world", { mass = 0.0 })

  for i = 1, 3 do
    local p_name = (i == 1) and "world" or ("link" .. (i - 1))
    local origin_z = (i == 1) and 0.0 or -l[i - 1]
    local mi = m[i]
    local li = l[i]
    local I = (1.0 / 12.0) * mi * li * li

    b:add_link({
      name = "link" .. i,
      parent = p_name,
      joint_name = "joint" .. i,
      joint_type = "revolute",
      axis = vec3.unit_y(),
      origin = { xyz = vec3.new(0.0, 0.0, origin_z), rpy = vec3.zero() },
      mass = mi,
      com = vec3.new(0.0, 0.0, -0.5 * li),
      inertia = { ixx = I, iyy = I, izz = 0.001 }
    })
  end

  local model, meta = b:build()
  meta.masses = m
  meta.lengths = l
  return model, meta
end

return pendulums
