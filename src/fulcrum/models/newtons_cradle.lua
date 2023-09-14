local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")
local shapes = require("fulcrum.geometry.shapes")

local newtons_cradle = {}

function newtons_cradle.create(params)
  params = params or {}
  local num_balls = params.num_balls or 5
  local ball_radius = params.radius or 0.05
  local ball_mass = params.mass or 0.2
  local string_length = params.length or 0.5
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("newtons_cradle", gravity)
  b:add_base("header_rail", { mass = 0.0 })

  local ball_shapes = {}
  local I_sphere = (2.0 / 5.0) * ball_mass * ball_radius * ball_radius

  for i = 1, num_balls do
    -- Balls aligned along X axis with spacing 2 * radius
    local x_pos = (i - (num_balls + 1) * 0.5) * (2.0 * ball_radius)

    b:add_link({
      name = "ball_" .. i,
      parent = "header_rail",
      joint_name = "pivot_" .. i,
      joint_type = "revolute",
      axis = vec3.unit_y(),
      origin = { xyz = vec3.new(x_pos, 0.0, 0.0), rpy = vec3.zero() },
      mass = ball_mass,
      com = vec3.new(0.0, 0.0, -string_length),
      inertia = { ixx = I_sphere, iyy = I_sphere, izz = I_sphere }
    })

    table.insert(ball_shapes, shapes.sphere(ball_radius, vec3.new(x_pos, 0.0, -string_length)))
  end

  local model, meta = b:build()
  meta.num_balls = num_balls
  meta.radius = ball_radius
  meta.mass = ball_mass
  meta.length = string_length
  meta.shapes = ball_shapes
  return model, meta
end

return newtons_cradle
