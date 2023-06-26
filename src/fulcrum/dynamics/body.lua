local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")

local body = {}
local body_mt = { __index = body }

function body.new(params)
  params = params or {}
  local mass = params.mass or 1.0
  local com = params.com or vec3.zero()
  local I_rot = params.I_rot or mat3.from_diag(1.0, 1.0, 1.0)

  local I_sp = spatial.inertia(mass, com, I_rot)

  local b = {
    id = params.id or 1,
    name = params.name or "body",
    parent_id = params.parent_id or 0,
    joint = params.joint,
    X_T = params.X_T or matrix.eye(6),
    mass = mass,
    com = com,
    I_rot = I_rot,
    I_spatial = I_sp
  }

  return setmetatable(b, body_mt)
end

return body
