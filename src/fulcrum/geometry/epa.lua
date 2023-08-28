local vec3 = require("fulcrum.math.vec3")
local gjk = require("fulcrum.geometry.gjk")

local epa = {}

local function get_face_normal_and_dist(a, b, c)
  local ab = b - a
  local ac = c - a
  local n = vec3.cross(ab, ac)
  local norm_n = vec3.norm(n)

  if norm_n < 1e-12 then
    return vec3.unit_z(), 0.0
  end

  local unit_n = n * (1.0 / norm_n)
  local d = vec3.dot(unit_n, a)

  if d < 0 then
    -- Ensure normal points away from origin
    unit_n = -unit_n
    d = -d
  end

  return unit_n, d
end

function epa.penetration(shape_a, shape_b, simplex, max_iter, tol)
  -- Expanding Polytope Algorithm: computes penetration depth and normal
  max_iter = max_iter or 64
  tol = tol or 1e-4

  -- Ensure simplex is a full 4-vertex tetrahedron enclosing origin
  if not simplex or #simplex < 4 then
    local hit, g_simplex = gjk.intersect(shape_a, shape_b)
    if not hit then
      return false, 0.0, vec3.zero()
    end
    simplex = g_simplex
  end

  -- Initial polytope faces from tetrahedron
  -- Faces array stores: { a, b, c, normal, dist }
  local faces = {}
  local function add_face(i1, i2, i3)
    local a, b, c = simplex[i1], simplex[i2], simplex[i3]
    local n, d = get_face_normal_and_dist(a, b, c)
    table.insert(faces, {
      v1 = a, v2 = b, v3 = c,
      normal = n,
      dist = d
    })
  end

  add_face(1, 2, 3)
  add_face(1, 3, 4)
  add_face(1, 4, 2)
  add_face(2, 4, 3)

  local best_face = nil
  for _ = 1, max_iter do
    -- Find face closest to origin
    local min_dist = 1e30
    local min_idx = 1

    for i = 1, #faces do
      if faces[i].dist < min_dist then
        min_dist = faces[i].dist
        min_idx = i
      end
    end

    best_face = faces[min_idx]

    -- Query support point in direction of closest face normal
    local sa = shape_a.support(best_face.normal)
    local sb = shape_b.support(-best_face.normal)
    local p = sa - sb

    local p_dist = vec3.dot(p, best_face.normal)
    if (p_dist - best_face.dist) < tol then
      -- Convergence: boundary reached
      return true, best_face.dist, best_face.normal, sa, sb
    end

    -- Find horizon edges of visible faces from p
    local edges = {}
    local function add_edge(ea, eb)
      for i = 1, #edges do
        local e = edges[i]
        if (vec3.dist_sq(e.a, eb) < 1e-12 and vec3.dist_sq(e.b, ea) < 1e-12) or
           (vec3.dist_sq(e.a, ea) < 1e-12 and vec3.dist_sq(e.b, eb) < 1e-12) then
          -- Shared edge: remove
          table.remove(edges, i)
          return
        end
      end
      table.insert(edges, { a = ea, b = eb })
    end

    local i = 1
    while i <= #faces do
      local f = faces[i]
      if vec3.dot(f.normal, p - f.v1) > 0 then
        -- Face is visible from new point p
        add_edge(f.v1, f.v2)
        add_edge(f.v2, f.v3)
        add_edge(f.v3, f.v1)
        table.remove(faces, i)
      else
        i = i + 1
      end
    end

    -- Construct new faces connecting horizon edges to p
    for _, edge in ipairs(edges) do
      local n, d = get_face_normal_and_dist(edge.a, edge.b, p)
      table.insert(faces, {
        v1 = edge.a, v2 = edge.b, v3 = p,
        normal = n,
        dist = d
      })
    end
  end

  if best_face then
    return true, best_face.dist, best_face.normal
  else
    return false, 0.0, vec3.zero()
  end
end

return epa
