local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local builder = require("fulcrum.models.builder")

local dzhanibekov = {}

function dzhanibekov.create(params)
  params = params or {}
  local mass = params.mass or 2.0
  -- Distinct principal moments: I1 < I2 < I3 (intermediate axis is Y)
  local I1 = params.I1 or 1.0
  local I2 = params.I2 or 2.0
  local I3 = params.I3 or 3.0

  local b = builder.new("dzhanibekov_top", vec3.zero())
  b:add_base("space", { mass = 0.0 })

  b:add_link({
    name = "wingnut",
    parent = "space",
    joint_name = "free_floating",
    joint_type = "floating_base",
    origin = { xyz = vec3.zero(), rpy = vec3.zero() },
    mass = mass,
    com = vec3.zero(),
    inertia = { ixx = I1, iyy = I2, izz = I3 }
  })

  local model, meta = b:build()
  meta.mass = mass
  meta.I1, meta.I2, meta.I3 = I1, I2, I3
  meta.I_diag = vec3.new(I1, I2, I3)
  return model, meta
end

return dzhanibekov
