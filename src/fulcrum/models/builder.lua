local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")
local joint_mod = require("fulcrum.dynamics.joint")
local model_mod = require("fulcrum.dynamics.model")

local builder = {}
local builder_mt = { __index = builder }

function builder.new(name, gravity)
  local b = {
    name = name or "robot",
    gravity = gravity or vec3.new(0.0, 0.0, -9.81),
    links = {},
    link_order = {},
    root_link = nil
  }
  return setmetatable(b, builder_mt)
end

local function parse_vec3(v, default)
  if not v then return default or vec3.zero() end
  if type(v) == "table" and v.x and v.y and v.z then
    return vec3.new(v.x, v.y, v.z)
  elseif type(v) == "table" and v[1] and v[2] and v[3] then
    return vec3.new(v[1], v[2], v[3])
  end
  return default or vec3.zero()
end

function builder:add_base(name, params)
  params = params or {}
  name = name or "base_link"
  self.root_link = name

  self.links[name] = {
    name = name,
    parent = nil,
    is_base = true,
    mass = params.mass or 0.0,
    com = parse_vec3(params.com),
    I_rot = params.I_rot or mat3.zero(),
    visual = params.visual,
    collision = params.collision
  }
  table.insert(self.link_order, name)
  return self
end

function builder:add_link(params)
  local name = params.name
  if not name then error("Link name is required") end

  local j_type = params.joint_type or "revolute"
  local j_axis = parse_vec3(params.axis, vec3.unit_z())
  local origin_xyz = parse_vec3(params.origin and params.origin.xyz or params.xyz)
  local origin_rpy = parse_vec3(params.origin and params.origin.rpy or params.rpy)

  local mass = params.mass or 1.0
  local com = parse_vec3(params.com)

  local I_rot
  if params.inertia then
    local I = params.inertia
    if type(I) == "table" and I[1] and type(I[1]) == "table" then
      I_rot = I
    else
      local ixx = I.ixx or I.Ixx or 0.1
      local iyy = I.iyy or I.Iyy or 0.1
      local izz = I.izz or I.Izz or 0.1
      local ixy = I.ixy or I.Ixy or 0.0
      local ixz = I.ixz or I.Ixz or 0.0
      local iyz = I.iyz or I.Iyz or 0.0
      I_rot = mat3.from_values(
        ixx, ixy, ixz,
        ixy, iyy, iyz,
        ixz, iyz, izz
      )
    end
  else
    I_rot = mat3.from_diag(0.1, 0.1, 0.1)
  end

  self.links[name] = {
    name = name,
    parent = params.parent,
    joint_name = params.joint_name or (name .. "_joint"),
    joint_type = j_type,
    axis = j_axis,
    origin_xyz = origin_xyz,
    origin_rpy = origin_rpy,
    mass = mass,
    com = com,
    I_rot = I_rot,
    limits = params.limits,
    visual = params.visual,
    collision = params.collision
  }
  table.insert(self.link_order, name)
  return self
end

function builder:build()
  local mdl = model_mod.new({ gravity = self.gravity })

  local link_id_map = {}
  local body_count = 0

  -- Topological sorting: start from base, then add children whose parents exist
  local added = {}
  local queue = {}

  if self.root_link and self.links[self.root_link] then
    local root = self.links[self.root_link]
    if root.mass > 0.0 then
      body_count = body_count + 1
      link_id_map[root.name] = body_count
      local j = joint_mod.fixed(root.name .. "_joint")
      mdl:add_body({
        name = root.name,
        parent_id = 0,
        joint = j,
        X_T = matrix.eye(6),
        mass = root.mass,
        com = root.com,
        I_rot = root.I_rot
      })
    else
      link_id_map[root.name] = 0
    end
    added[root.name] = true
  end

  -- Repeated pass to resolve parent dependencies
  local remaining = #self.link_order
  while remaining > 0 do
    local progress = false

    for _, name in ipairs(self.link_order) do
      if not added[name] then
        local lk = self.links[name]
        local p_name = lk.parent

        if not p_name or added[p_name] then
          local p_id = p_name and link_id_map[p_name] or 0

          -- Compute fixed coordinate transform X_T = spatial.transform(R, p)
          local R = mat3.from_euler_xyz(lk.origin_rpy.x, lk.origin_rpy.y, lk.origin_rpy.z)
          local X_T = spatial.transform(R, lk.origin_xyz)

          -- Instantiate joint
          local j
          local jt = lk.joint_type
          if jt == "revolute" or jt == "continuous" then
            j = joint_mod.revolute(lk.axis, lk.joint_name)
          elseif jt == "prismatic" then
            j = joint_mod.prismatic(lk.axis, lk.joint_name)
          elseif jt == "spherical" then
            j = joint_mod.spherical(lk.joint_name)
          elseif jt == "floating_base" then
            j = joint_mod.floating_base(lk.joint_name)
          else
            j = joint_mod.fixed(lk.joint_name)
          end

          if lk.limits then
            if lk.limits.lower then j.limit_lower = lk.limits.lower end
            if lk.limits.upper then j.limit_upper = lk.limits.upper end
            if lk.limits.effort then j.limit_effort = lk.limits.effort end
            if lk.limits.velocity then j.limit_velocity = lk.limits.velocity end
          end

          body_count = body_count + 1
          link_id_map[name] = body_count

          mdl:add_body({
            name = name,
            parent_id = p_id,
            joint = j,
            X_T = X_T,
            mass = lk.mass,
            com = lk.com,
            I_rot = lk.I_rot
          })

          added[name] = true
          progress = true
        end
      end
    end

    if not progress then
      error("Cyclic or disconnected link structure detected in robot model")
    end

    local count = 0
    for _, name in ipairs(self.link_order) do
      if not added[name] then count = count + 1 end
    end
    remaining = count
  end

  return mdl, {
    link_to_id = link_id_map,
    name = self.name
  }
end

return builder
