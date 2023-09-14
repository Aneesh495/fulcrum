local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local builder = require("fulcrum.models.builder")

local cartpole = {}

function cartpole.create(params)
  params = params or {}
  local cart_mass = params.cart_mass or 1.0
  local pole_mass = params.pole_mass or 0.1
  local pole_length = params.pole_length or 1.0
  local gravity = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local b = builder.new("cartpole", gravity)
  b:add_base("world", { mass = 0.0 })

  -- Cart: prismatic joint along X axis
  b:add_link({
    name = "cart",
    parent = "world",
    joint_name = "cart_slider",
    joint_type = "prismatic",
    axis = vec3.unit_x(),
    origin = { xyz = vec3.zero(), rpy = vec3.zero() },
    mass = cart_mass,
    com = vec3.zero(),
    inertia = { ixx = 0.1, iyy = 0.1, izz = 0.1 }
  })

  -- Pole: revolute joint rotating around Y axis
  local half_l = pole_length * 0.5
  local I_pole = (1.0 / 12.0) * pole_mass * pole_length * pole_length
  b:add_link({
    name = "pole",
    parent = "cart",
    joint_name = "pole_hinge",
    joint_type = "revolute",
    axis = vec3.unit_y(),
    origin = { xyz = vec3.zero(), rpy = vec3.zero() },
    mass = pole_mass,
    com = vec3.new(0.0, 0.0, half_l),
    inertia = { ixx = I_pole, iyy = I_pole, izz = 0.001 }
  })

  local model, meta = b:build()
  meta.cart_mass = cart_mass
  meta.pole_mass = pole_mass
  meta.pole_length = pole_length
  return model, meta
end

return cartpole
