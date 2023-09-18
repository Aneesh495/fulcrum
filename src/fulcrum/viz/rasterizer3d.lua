local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local terminal = require("fulcrum.viz.terminal")

local rasterizer3d = {}
local rast_mt = { __index = rasterizer3d }

function rasterizer3d.new(width, height)
  width = width or 80
  height = height or 32

  local r = {
    width = width,
    height = height,
    aspect = (width / height) * 0.5, -- Adjust for non-square terminal character cells
    fov = math.rad(60.0),
    near = 0.1,
    far = 100.0,
    cam_pos = vec3.new(3.0, -3.0, 2.0),
    cam_target = vec3.new(0.0, 0.0, 0.0),
    cam_up = vec3.unit_z(),
    char_buf = {},
    color_buf = {},
    depth_buf = {}
  }

  setmetatable(r, rast_mt)
  r:clear()
  return r
end

function rasterizer3d:clear()
  local w, h = self.width, self.height
  for y = 1, h do
    self.char_buf[y] = {}
    self.color_buf[y] = {}
    self.depth_buf[y] = {}
    for x = 1, w do
      self.char_buf[y][x] = " "
      self.color_buf[y][x] = nil
      self.depth_buf[y][x] = -1e30
    end
  end
end

function rasterizer3d:set_camera(pos, target, up)
  self.cam_pos = pos or self.cam_pos
  self.cam_target = target or self.cam_target
  self.cam_up = up or self.cam_up
end

function rasterizer3d:project(p_world)
  -- Transform from world to camera coordinate frame
  local f = vec3.unit(self.cam_target - self.cam_pos) -- forward
  local r = vec3.unit(vec3.cross(f, self.cam_up))     -- right
  local u = vec3.cross(r, f)                          -- true up

  local rel = p_world - self.cam_pos
  local xc = vec3.dot(rel, r)
  local yc = vec3.dot(rel, u)
  local zc = vec3.dot(rel, f) -- depth along camera optical axis

  if zc <= self.near then
    return nil, nil, nil
  end

  local tan_half = math.tan(self.fov * 0.5)
  local xs_norm = xc / (zc * tan_half * self.aspect)
  local ys_norm = yc / (zc * tan_half)

  local sx = math.floor((xs_norm + 1.0) * 0.5 * self.width) + 1
  local sy = math.floor((1.0 - (ys_norm + 1.0) * 0.5) * self.height) + 1

  return sx, sy, 1.0 / zc
end

function rasterizer3d:set_pixel(x, y, inv_z, char, color)
  if x < 1 or x > self.width or y < 1 or y > self.height then
    return
  end

  if inv_z >= self.depth_buf[y][x] then
    self.depth_buf[y][x] = inv_z
    self.char_buf[y][x] = char or "*"
    self.color_buf[y][x] = color
  end
end

function rasterizer3d:draw_line_3d(p1, p2, char, color)
  local sx1, sy1, iz1 = self:project(p1)
  local sx2, sy2, iz2 = self:project(p2)

  if not sx1 or not sx2 then return end

  local dx = math.abs(sx2 - sx1)
  local dy = math.abs(sy2 - sy1)
  local sx = (sx1 < sx2) and 1 or -1
  local sy = (sy1 < sy2) and 1 or -1
  local err = dx - dy

  local x, y = sx1, sy1
  local total_steps = math.max(dx, dy)
  local step = 0

  while true do
    local alpha = (total_steps > 0) and (step / total_steps) or 0.0
    local iz = (1.0 - alpha) * iz1 + alpha * iz2
    self:set_pixel(x, y, iz, char, color)

    if x == sx2 and y == sy2 then break end

    local e2 = 2 * err
    if e2 > -dy then
      err = err - dy
      x = x + sx
    end
    if e2 < dx then
      err = err + dx
      y = y + sy
    end
    step = step + 1
  end
end

function rasterizer3d:draw_axes(origin, R, scale)
  scale = scale or 0.5
  origin = origin or vec3.zero()
  R = R or mat3.identity()

  local ax_x = origin + (R * vec3.new(scale, 0, 0))
  local ax_y = origin + (R * vec3.new(0, scale, 0))
  local ax_z = origin + (R * vec3.new(0, 0, scale))

  self:draw_line_3d(origin, ax_x, "X", terminal.rgb_fg(255, 60, 60))
  self:draw_line_3d(origin, ax_y, "Y", terminal.rgb_fg(60, 255, 60))
  self:draw_line_3d(origin, ax_z, "Z", terminal.rgb_fg(60, 160, 255))
end

function rasterizer3d:draw_grid(size, step)
  size = size or 2.0
  step = step or 0.5
  local color = terminal.rgb_fg(80, 80, 80)

  local x = -size
  while x <= size + 1e-6 do
    self:draw_line_3d(vec3.new(x, -size, 0), vec3.new(x, size, 0), "·", color)
    x = x + step
  end

  local y = -size
  while y <= size + 1e-6 do
    self:draw_line_3d(vec3.new(-size, y, 0), vec3.new(size, y, 0), "·", color)
    y = y + step
  end
end

function rasterizer3d:draw_box_wireframe(center, half_extents, R, color)
  center = center or vec3.zero()
  R = R or mat3.identity()
  local h = half_extents or vec3.new(0.5, 0.5, 0.5)

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

  for i = 1, 8 do
    corners[i] = center + (R * corners[i])
  end

  local edges = {
    {1, 2}, {2, 3}, {3, 4}, {4, 1},
    {5, 6}, {6, 7}, {7, 8}, {8, 5},
    {1, 5}, {2, 6}, {3, 7}, {4, 8}
  }

  for _, e in ipairs(edges) do
    self:draw_line_3d(corners[e[1]], corners[e[2]], "#", color)
  end
end

function rasterizer3d:render_to_string()
  local lines = {}
  for y = 1, self.height do
    local line_parts = {}
    local cur_color = nil

    for x = 1, self.width do
      local col = self.color_buf[y][x]
      local ch = self.char_buf[y][x]

      if col ~= cur_color then
        if col then
          table.insert(line_parts, col)
        else
          table.insert(line_parts, terminal.RESET)
        end
        cur_color = col
      end
      table.insert(line_parts, ch)
    end

    if cur_color then
      table.insert(line_parts, terminal.RESET)
    end
    table.insert(lines, table.concat(line_parts))
  end

  return table.concat(lines, "\n")
end

return rasterizer3d
