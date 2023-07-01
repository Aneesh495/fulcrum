local vec3 = require("fulcrum.math.vec3")
local matrix = require("fulcrum.math.matrix")
local body = require("fulcrum.dynamics.body")

local model = {}
local model_mt = { __index = model }

function model.new(params)
  params = params or {}
  local g = params.gravity or vec3.new(0.0, 0.0, -9.81)

  local m = {
    bodies = {},
    num_bodies = 0,
    dof = 0,
    gravity = { 0.0, 0.0, 0.0, g.x, g.y, g.z }, -- spatial acceleration of gravity
    q_indices = {},
    v_indices = {}
  }

  return setmetatable(m, model_mt)
end

function model:add_body(params)
  self.num_bodies = self.num_bodies + 1
  params.id = self.num_bodies

  local b = body.new(params)
  table.insert(self.bodies, b)

  local dof = b.joint and b.joint.dof or 0
  local q_start = self.dof + 1
  local q_end = self.dof + dof
  self.q_indices[b.id] = { start = q_start, finish = q_end, dof = dof }
  self.v_indices[b.id] = { start = q_start, finish = q_end, dof = dof }
  self.dof = self.dof + dof

  return b.id
end

function model:get_body(id)
  return self.bodies[id]
end

return model
