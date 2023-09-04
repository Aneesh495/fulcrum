local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")

local shapes = {}

function shapes.sphere(radius, center)
  local r = radius or 0.5
  local c = center or vec3.zero()
  return {
    type = "sphere",
    radius = r,
    center = c,
    support = function(d)
      local u = vec3.unit(d)
      return c + (u * r)
    end,
    aabb = function()
      local ext = vec3.new(r, r, r)
      return c - ext, c + ext
    end
  }
end

function shapes.box(half_extents, transform)
  local h = half_extents or vec3.new(0.5, 0.5, 0.5)
  local T = transform -- { R, p }
  return {
    type = "box",
    half_extents = h,
    transform = T,
    support = function(d)
      local d_local = d
      if T and T.R then
        local Rt = mat3.transpose(T.R)
        d_local = Rt * d
      end

      local sx = (d_local.x >= 0) and h.x or -h.x
      local sy = (d_local.y >= 0) and h.y or -h.y
      local sz = (d_local.z >= 0) and h.z or -h.z
      local v_local = vec3.new(sx, sy, sz)

      if T then
        local v_world = (T.R and (T.R * v_local) or v_local) + (T.p or vec3.zero())
        return v_world
      else
        return v_local
      end
    end,
    aabb = function()
      local p = (T and T.p) or vec3.zero()
      local R = (T and T.R) or mat3.identity()
      -- Ex = |R11|*hx + |R12|*hy + |R13|*hz, etc.
      local ex = math.abs(R[1][1])*h.x + math.abs(R[1][2])*h.y + math.abs(R[1][3])*h.z
      local ey = math.abs(R[2][1])*h.x + math.abs(R[2][2])*h.y + math.abs(R[2][3])*h.z
      local ez = math.abs(R[3][1])*h.x + math.abs(R[3][2])*h.y + math.abs(R[3][3])*h.z
      local ext = vec3.new(ex, ey, ez)
      return p - ext, p + ext
    end
  }
end

function shapes.capsule(radius, p1, p2)
  local r = radius or 0.2
  local a = p1 or vec3.new(0, 0, -0.5)
  local b = p2 or vec3.new(0, 0, 0.5)
  return {
    type = "capsule",
    radius = r,
    p1 = a,
    p2 = b,
    support = function(d)
      local u = vec3.unit(d)
      local da = vec3.dot(a, d)
      local db = vec3.dot(b, d)
      local base = (da >= db) and a or b
      return base + (u * r)
    end,
    aabb = function()
      local min_x = math.min(a.x, b.x) - r
      local min_y = math.min(a.y, b.y) - r
      local min_z = math.min(a.z, b.z) - r
      local max_x = math.max(a.x, b.x) + r
      local max_y = math.max(a.y, b.y) + r
      local max_z = math.max(a.z, b.z) + r
      return vec3.new(min_x, min_y, min_z), vec3.new(max_x, max_y, max_z)
    end
  }
end

function shapes.cylinder(radius, half_height, transform)
  local r = radius or 0.3
  local h = half_height or 0.5
  local T = transform
  return {
    type = "cylinder",
    radius = r,
    half_height = h,
    transform = T,
    support = function(d)
      local d_local = d
      if T and T.R then
        local Rt = mat3.transpose(T.R)
        d_local = Rt * d
      end

      local sigma = math.sqrt(d_local.x * d_local.x + d_local.y * d_local.y)
      local px, py
      if sigma > 1e-12 then
        px = r * d_local.x / sigma
        py = r * d_local.y / sigma
      else
        px = 0.0
        py = 0.0
      end
      local pz = (d_local.z >= 0) and h or -h
      local v_local = vec3.new(px, py, pz)

      if T then
        return (T.R and (T.R * v_local) or v_local) + (T.p or vec3.zero())
      else
        return v_local
      end
    end
  }
end

function shapes.ellipsoid(semi_axes, transform)
  local a, b, c = semi_axes.x, semi_axes.y, semi_axes.z
  local T = transform
  return {
    type = "ellipsoid",
    semi_axes = semi_axes,
    support = function(d)
      local d_local = d
      if T and T.R then
        local Rt = mat3.transpose(T.R)
        d_local = Rt * d
      end

      -- S(d) = [ a^2 dx, b^2 dy, c^2 dz ] / sqrt(a^2 dx^2 + b^2 dy^2 + c^2 dz^2)
      local a2 = a * a
      local b2 = b * b
      local c2 = c * c

      local kx = a2 * d_local.x
      local ky = b2 * d_local.y
      local kz = c2 * d_local.z

      local denom = math.sqrt(d_local.x * kx + d_local.y * ky + d_local.z * kz)
      local v_local
      if denom > 1e-15 then
        local inv = 1.0 / denom
        v_local = vec3.new(kx * inv, ky * inv, kz * inv)
      else
        v_local = vec3.zero()
      end

      if T then
        return (T.R and (T.R * v_local) or v_local) + (T.p or vec3.zero())
      else
        return v_local
      end
    end
  }
end

function shapes.polyhedron(vertices, transform)
  local verts = vertices -- list of vec3
  local T = transform
  return {
    type = "polyhedron",
    vertices = verts,
    support = function(d)
      local max_dot = -1e30
      local best_v = verts[1] or vec3.zero()

      for i = 1, #verts do
        local v = verts[i]
        local v_world = T and ((T.R and (T.R * v) or v) + (T.p or vec3.zero())) or v
        local dot = vec3.dot(v_world, d)
        if dot > max_dot then
          max_dot = dot
          best_v = v_world
        end
      end

      return best_v
    end
  }
end

function shapes.translate(shape, delta)
  delta = delta or vec3.zero()
  local t = shape.type
  if t == "sphere" then
    return shapes.sphere(shape.radius, shape.center + delta)
  elseif t == "box" then
    local p = (shape.transform and shape.transform.p) or vec3.zero()
    local R = (shape.transform and shape.transform.R) or mat3.identity()
    return shapes.box(shape.half_extents, { R = R, p = p + delta })
  elseif t == "capsule" then
    local p = (shape.transform and shape.transform.p) or vec3.zero()
    local R = (shape.transform and shape.transform.R) or mat3.identity()
    return shapes.capsule(shape.radius, shape.half_height, { R = R, p = p + delta })
  elseif t == "cylinder" then
    local p = (shape.transform and shape.transform.p) or vec3.zero()
    local R = (shape.transform and shape.transform.R) or mat3.identity()
    return shapes.cylinder(shape.radius, shape.half_height, { R = R, p = p + delta })
  elseif t == "ellipsoid" then
    local p = (shape.transform and shape.transform.p) or vec3.zero()
    local R = (shape.transform and shape.transform.R) or mat3.identity()
    return shapes.ellipsoid(shape.radii, { R = R, p = p + delta })
  elseif t == "polyhedron" then
    local p = (shape.transform and shape.transform.p) or vec3.zero()
    local R = (shape.transform and shape.transform.R) or mat3.identity()
    return shapes.polyhedron(shape.vertices, { R = R, p = p + delta })
  else
    return {
      type = "translated_" .. tostring(t),
      support = function(d)
        return shape.support(d) + delta
      end,
      translated = function(self, d2)
        return shapes.translate(self, d2)
      end
    }
  end
end

return shapes
