local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")

local obj_exporter = {}

function obj_exporter.box_mesh(half_extents, transform)
  local h = half_extents or vec3.new(0.5, 0.5, 0.5)
  local T = transform or {}
  local p = T.p or vec3.zero()
  local R = T.R or mat3.identity()

  local corners = {
    vec3.new(-h.x, -h.y, -h.z),
    vec3.new( h.x, -h.y, -h.z),
    vec3.new( h.x,  h.y, -h.z),
    vec3.new(-h.x,  h.y, -h.z),
    vec3.new(-h.x, -h.y,  h.z),
    vec3.new( h.x, -h.y,  h.z),
    vec3.new( h.x,  h.y,  h.z),
    vec3.new(-h.x,  h.y,  h.z)
  }

  local verts = {}
  for i = 1, 8 do
    verts[i] = p + (R * corners[i])
  end

  local faces = {
    {1, 2, 3}, {1, 3, 4}, -- bottom
    {5, 7, 6}, {5, 8, 7}, -- top
    {1, 6, 2}, {1, 5, 6}, -- front
    {2, 7, 3}, {2, 6, 7}, -- right
    {3, 8, 4}, {3, 7, 8}, -- back
    {4, 5, 1}, {4, 8, 5}  -- left
  }

  return verts, faces
end

function obj_exporter.to_obj_string(meshes)
  local lines = {
    "# Fulcrum Wavefront OBJ Exporter",
    "# Units: meters, right-handed Z-up coordinate system"
  }

  local v_offset = 0

  for idx, mesh in ipairs(meshes) do
    local name = mesh.name or ("object_" .. idx)
    table.insert(lines, "o " .. name)

    for _, v in ipairs(mesh.vertices) do
      table.insert(lines, string.format("v %.6f %.6f %.6f", v.x, v.y, v.z))
    end

    for _, f in ipairs(mesh.faces) do
      table.insert(lines, string.format("f %d %d %d",
        f[1] + v_offset,
        f[2] + v_offset,
        f[3] + v_offset
      ))
    end

    v_offset = v_offset + #mesh.vertices
  end

  return table.concat(lines, "\n")
end

function obj_exporter.write_file(filepath, meshes)
  local str = obj_exporter.to_obj_string(meshes)
  local f, err = io.open(filepath, "w")
  if not f then
    error(string.format("Failed to open OBJ file '%s': %s", tostring(filepath), tostring(err)))
  end
  f:write(str)
  f:close()
  return true
end

return obj_exporter
