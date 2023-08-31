local vec3 = require("fulcrum.math.vec3")

local bvh = {}

-- AABB representation and utilities
local aabb = {}
aabb.__index = aabb

function aabb.new(min_pt, max_pt)
  return setmetatable({
    min = min_pt or vec3.new(0, 0, 0),
    max = max_pt or vec3.new(0, 0, 0)
  }, aabb)
end

function aabb.from_point(pt, margin)
  margin = margin or 0.0
  local m = vec3.new(margin, margin, margin)
  return aabb.new(pt - m, pt + m)
end

function aabb.surface_area(box)
  local d = box.max - box.min
  local dx = math.max(0.0, d.x)
  local dy = math.max(0.0, d.y)
  local dz = math.max(0.0, d.z)
  return 2.0 * (dx * dy + dy * dz + dz * dx)
end

function aabb.union(b1, b2)
  return aabb.new(
    vec3.new(
      math.min(b1.min.x, b2.min.x),
      math.min(b1.min.y, b2.min.y),
      math.min(b1.min.z, b2.min.z)
    ),
    vec3.new(
      math.max(b1.max.x, b2.max.x),
      math.max(b1.max.y, b2.max.y),
      math.max(b1.max.z, b2.max.z)
    )
  )
end

function aabb.intersects(b1, b2)
  if b1.max.x < b2.min.x or b1.min.x > b2.max.x then return false end
  if b1.max.y < b2.min.y or b1.min.y > b2.max.y then return false end
  if b1.max.z < b2.min.z or b1.min.z > b2.max.z then return false end
  return true
end

function aabb.contains_point(box, pt)
  return pt.x >= box.min.x and pt.x <= box.max.x and
         pt.y >= box.min.y and pt.y <= box.max.y and
         pt.z >= box.min.z and pt.z <= box.max.z
end

function aabb.ray_intersect(box, origin, dir, max_t)
  max_t = max_t or 1e30
  local tmin = -1e30
  local tmax = max_t

  local dims = { "x", "y", "z" }
  for _, d in ipairs(dims) do
    local inv_d = 1.0 / (math.abs(dir[d]) > 1e-12 and dir[d] or (dir[d] >= 0 and 1e-12 or -1e-12))
    local t0 = (box.min[d] - origin[d]) * inv_d
    local t1 = (box.max[d] - origin[d]) * inv_d

    if inv_d < 0.0 then
      t0, t1 = t1, t0
    end

    tmin = math.max(tmin, t0)
    tmax = math.min(tmax, t1)

    if tmax < tmin or tmax < 0.0 then
      return false, 0.0
    end
  end

  return true, math.max(0.0, tmin), tmax
end

bvh.aabb = aabb

-- Dynamic BVH Tree with Surface Area Heuristic
local tree = {}
tree.__index = tree

function bvh.new_tree()
  return setmetatable({
    root = nil,
    node_count = 0,
    leaf_nodes = {}
  }, tree)
end

local function new_node(box, user_id, user_data)
  return {
    box = box,
    user_id = user_id,
    user_data = user_data,
    left = nil,
    right = nil,
    parent = nil,
    is_leaf = (user_id ~= nil)
  }
end

function tree:insert(user_id, box, user_data)
  local leaf = new_node(box, user_id, user_data)
  self.leaf_nodes[user_id] = leaf
  self.node_count = self.node_count + 1

  if not self.root then
    self.root = leaf
    return leaf
  end

  -- Stage 1: find best sibling for leaf using Surface Area Heuristic
  local current = self.root
  while not current.is_leaf do
    local combined_box = aabb.union(current.box, box)
    local combined_cost = 2.0 * aabb.surface_area(combined_box)

    -- Minimum cost of pushing leaf down tree
    local inherited_cost = 2.0 * (aabb.surface_area(combined_box) - aabb.surface_area(current.box))

    -- Cost of descending left
    local cost_left
    local left_box = aabb.union(current.left.box, box)
    if current.left.is_leaf then
      cost_left = aabb.surface_area(left_box) + inherited_cost
    else
      local old_area = aabb.surface_area(current.left.box)
      local new_area = aabb.surface_area(left_box)
      cost_left = (new_area - old_area) + inherited_cost
    end

    -- Cost of descending right
    local cost_right
    local right_box = aabb.union(current.right.box, box)
    if current.right.is_leaf then
      cost_right = aabb.surface_area(right_box) + inherited_cost
    else
      local old_area = aabb.surface_area(current.right.box)
      local new_area = aabb.surface_area(right_box)
      cost_right = (new_area - old_area) + inherited_cost
    end

    -- Descend according to minimum cost
    if combined_cost < cost_left and combined_cost < cost_right then
      break
    end

    if cost_left < cost_right then
      current = current.left
    else
      current = current.right
    end
  end

  local sibling = current

  -- Stage 2: create new parent node
  local old_parent = sibling.parent
  local new_parent = new_node(aabb.union(sibling.box, box), nil, nil)
  new_parent.parent = old_parent
  new_parent.left = sibling
  new_parent.right = leaf

  sibling.parent = new_parent
  leaf.parent = new_parent

  if not old_parent then
    self.root = new_parent
  else
    if old_parent.left == sibling then
      old_parent.left = new_parent
    else
      old_parent.right = new_parent
    end
  end

  -- Stage 3: refit bounding boxes up ancestor chain
  local ancestor = new_parent
  while ancestor do
    ancestor.box = aabb.union(ancestor.left.box, ancestor.right.box)
    ancestor = ancestor.parent
  end

  return leaf
end

function tree:remove(user_id)
  local leaf = self.leaf_nodes[user_id]
  if not leaf then return false end

  self.leaf_nodes[user_id] = nil
  self.node_count = self.node_count - 1

  if leaf == self.root then
    self.root = nil
    return true
  end

  local parent = leaf.parent
  local grandparent = parent.parent
  local sibling = (parent.left == leaf) and parent.right or parent.left

  if not grandparent then
    self.root = sibling
    sibling.parent = nil
  else
    if grandparent.left == parent then
      grandparent.left = sibling
    else
      grandparent.right = sibling
    end
    sibling.parent = grandparent

    -- Refit ancestors
    local ancestor = grandparent
    while ancestor do
      ancestor.box = aabb.union(ancestor.left.box, ancestor.right.box)
      ancestor = ancestor.parent
    end
  end

  return true
end

function tree:query_aabb(query_box, callback)
  if not self.root then return end

  local stack = { self.root }
  while #stack > 0 do
    local node = table.remove(stack)
    if aabb.intersects(node.box, query_box) then
      if node.is_leaf then
        callback(node.user_id, node.user_data, node.box)
      else
        if node.left then table.insert(stack, node.left) end
        if node.right then table.insert(stack, node.right) end
      end
    end
  end
end

function tree:query_pairs(callback)
  if not self.root or self.root.is_leaf then return end

  local function collide(n1, n2)
    if not aabb.intersects(n1.box, n2.box) then
      return
    end

    if n1.is_leaf and n2.is_leaf then
      callback(n1.user_id, n2.user_id, n1.user_data, n2.user_data)
      return
    end

    if n1.is_leaf then
      collide(n1, n2.left)
      collide(n1, n2.right)
    elseif n2.is_leaf then
      collide(n1.left, n2)
      collide(n1.right, n2)
    else
      -- Split larger node by surface area
      if aabb.surface_area(n1.box) > aabb.surface_area(n2.box) then
        collide(n1.left, n2)
        collide(n1.right, n2)
      else
        collide(n1, n2.left)
        collide(n1, n2.right)
      end
    end
  end

  local function self_collide(node)
    if node.is_leaf then return end
    collide(node.left, node.right)
    self_collide(node.left)
    self_collide(node.right)
  end

  self_collide(self.root)
end

function tree:raycast(origin, dir, max_dist)
  max_dist = max_dist or 1e30
  local hits = {}
  if not self.root then return hits end

  local stack = { self.root }
  while #stack > 0 do
    local node = table.remove(stack)
    local hit, tmin, tmax = aabb.ray_intersect(node.box, origin, dir, max_dist)
    if hit and tmin <= max_dist then
      if node.is_leaf then
        table.insert(hits, {
          user_id = node.user_id,
          user_data = node.user_data,
          tmin = tmin,
          tmax = tmax
        })
      else
        if node.left then table.insert(stack, node.left) end
        if node.right then table.insert(stack, node.right) end
      end
    end
  end

  table.sort(hits, function(a, b) return a.tmin < b.tmin end)
  return hits
end

return bvh
