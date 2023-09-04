local vec3 = require("fulcrum.math.vec3")

local gjk = {}

local function support_minkowski(shape_a, shape_b, dir)
  local sa = shape_a.support(dir)
  local sb = shape_b.support(-dir)
  return sa - sb, sa, sb
end

local function handle_simplex(simplex, pts_a, pts_b, dir)
  local n = #simplex

  if n == 2 then
    -- Line segment AB (A is newest point, B is older)
    local A = simplex[2]
    local B = simplex[1]
    local AB = B - A
    local AO = -A

    if vec3.dot(AB, AO) > 0 then
      -- Origin is in region of segment AB
      local new_dir = vec3.cross(vec3.cross(AB, AO), AB)
      if vec3.norm_sq(new_dir) < 1e-12 then
        -- Collinear: pick perpendicular
        new_dir = vec3.cross(AB, vec3.unit_x())
        if vec3.norm_sq(new_dir) < 1e-12 then
          new_dir = vec3.cross(AB, vec3.unit_y())
        end
      end
      dir:copy(new_dir)
      return false
    else
      -- Origin is past A
      simplex[1] = A
      pts_a[1] = pts_a[2]
      pts_b[1] = pts_b[2]
      table.remove(simplex, 2)
      table.remove(pts_a, 2)
      table.remove(pts_b, 2)
      dir:copy(AO)
      return false
    end

  elseif n == 3 then
    -- Triangle ABC (A is newest)
    local A = simplex[3]
    local B = simplex[2]
    local C = simplex[1]

    local AB = B - A
    local AC = C - A
    local AO = -A
    local ABC = vec3.cross(AB, AC)

    -- Check edge AB normal
    local AB_perp = vec3.cross(AB, ABC)
    if vec3.dot(AB_perp, AO) > 0 then
      -- Outside edge AB
      table.remove(simplex, 1) -- remove C
      table.remove(pts_a, 1)
      table.remove(pts_b, 1)
      dir:copy(vec3.cross(vec3.cross(AB, AO), AB))
      return false
    end

    -- Check edge AC normal
    local AC_perp = vec3.cross(ABC, AC)
    if vec3.dot(AC_perp, AO) > 0 then
      -- Outside edge AC
      table.remove(simplex, 2) -- remove B
      table.remove(pts_a, 2)
      table.remove(pts_b, 2)
      dir:copy(vec3.cross(vec3.cross(AC, AO), AC))
      return false
    end

    -- Origin is above or below triangle ABC
    if vec3.dot(ABC, AO) > 0 then
      dir:copy(ABC)
    else
      -- Flip winding
      simplex[1], simplex[2] = simplex[2], simplex[1]
      pts_a[1], pts_a[2] = pts_a[2], pts_a[1]
      pts_b[1], pts_b[2] = pts_b[2], pts_b[1]
      dir:copy(-ABC)
    end
    return false

  elseif n == 4 then
    -- Tetrahedron ABCD (A is newest)
    local A = simplex[4]
    local B = simplex[3]
    local C = simplex[2]
    local D = simplex[1]

    local AB = B - A
    local AC = C - A
    local AD = D - A
    local AO = -A

    local norm_ABC = vec3.cross(AB, AC)
    local norm_ACD = vec3.cross(AC, AD)
    local norm_ADB = vec3.cross(AD, AB)

    -- Check face ABC
    if vec3.dot(norm_ABC, AO) > 0 then
      table.remove(simplex, 1) -- remove D
      table.remove(pts_a, 1)
      table.remove(pts_b, 1)
      dir:copy(norm_ABC)
      return false
    end

    -- Check face ACD
    if vec3.dot(norm_ACD, AO) > 0 then
      table.remove(simplex, 3) -- remove B
      table.remove(pts_a, 3)
      table.remove(pts_b, 3)
      dir:copy(norm_ACD)
      return false
    end

    -- Check face ADB
    if vec3.dot(norm_ADB, AO) > 0 then
      table.remove(simplex, 2) -- remove C
      table.remove(pts_a, 2)
      table.remove(pts_b, 2)
      dir:copy(norm_ADB)
      return false
    end

    -- Origin is inside tetrahedron ABCD: collision detected!
    return true
  end

  return false
end

function gjk.intersect(shape_a, shape_b, max_iter)
  max_iter = max_iter or 40

  -- Initial search direction: difference between centers or arbitrary axis
  local d = vec3.new(1.0, 0.0, 0.0)

  local simplex = {}
  local pts_a = {}
  local pts_b = {}

  local first_pt, p_a, p_b = support_minkowski(shape_a, shape_b, d)
  table.insert(simplex, first_pt)
  table.insert(pts_a, p_a)
  table.insert(pts_b, p_b)

  d = -first_pt

  for _ = 1, max_iter do
    if vec3.norm_sq(d) < 1e-14 then
      -- Origin already on boundary
      return true, simplex, pts_a, pts_b
    end

    local new_pt, pa, pb = support_minkowski(shape_a, shape_b, d)
    if vec3.dot(new_pt, d) < 0 then
      -- Cannot pass origin: shapes are disjoint
      return false, simplex, pts_a, pts_b
    end

    table.insert(simplex, new_pt)
    table.insert(pts_a, pa)
    table.insert(pts_b, pb)

    if handle_simplex(simplex, pts_a, pts_b, d) then
      return true, simplex, pts_a, pts_b
    end
  end

  return false, simplex, pts_a, pts_b
end

local function closest_point_segment(a, b)
  local ab = b - a
  local d = vec3.dot(-a, ab)
  local len_sq = vec3.dot(ab, ab)
  if d <= 0 or len_sq < 1e-14 then
    return a, { 1, 0 }, { 1 }
  elseif d >= len_sq then
    return b, { 0, 1 }, { 2 }
  else
    local t = d / len_sq
    return a + ab * t, { 1 - t, t }, { 1, 2 }
  end
end

local function closest_point_triangle(a, b, c)
  local ab = b - a
  local ac = c - a
  local ap = -a

  local d1 = vec3.dot(ab, ap)
  local d2 = vec3.dot(ac, ap)
  if d1 <= 0 and d2 <= 0 then
    return a, { 1, 0, 0 }, { 1 }
  end

  local bp = -b
  local d3 = vec3.dot(ab, bp)
  local d4 = vec3.dot(ac, bp)
  if d3 >= 0 and d4 <= d3 then
    return b, { 0, 1, 0 }, { 2 }
  end

  local cp = -c
  local d5 = vec3.dot(ab, cp)
  local d6 = vec3.dot(ac, cp)
  if d6 >= 0 and d5 <= d6 then
    return c, { 0, 0, 1 }, { 3 }
  end

  local vc = d1 * d4 - d3 * d2
  if vc <= 0 and d1 >= 0 and d3 <= 0 then
    local v = d1 / (d1 - d3)
    return a + ab * v, { 1 - v, v, 0 }, { 1, 2 }
  end

  local vb = d5 * d2 - d1 * d6
  if vb <= 0 and d2 >= 0 and d6 <= 0 then
    local w = d2 / (d2 - d6)
    return a + ac * w, { 1 - w, 0, w }, { 1, 3 }
  end

  local va = d3 * d6 - d5 * d4
  if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0 then
    local w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
    return b + (c - b) * w, { 0, 1 - w, w }, { 2, 3 }
  end

  local denom = 1.0 / (va + vb + vc)
  local v = vb * denom
  local w = vc * denom
  local u = 1.0 - v - w
  return a * u + b * v + c * w, { u, v, w }, { 1, 2, 3 }
end

function gjk.distance(shape_a, shape_b, max_iter, tol)
  max_iter = max_iter or 30
  tol = tol or 1e-5

  local hit, simplex, pts_a, pts_b = gjk.intersect(shape_a, shape_b)
  if hit then
    return 0.0, nil, nil, vec3.zero()
  end

  -- Initial direction
  local d = vec3.unit_x()
  local pt, pa, pb = support_minkowski(shape_a, shape_b, d)
  local s = { pt }
  local sa = { pa }
  local sb = { pb }

  local last_dist = 1e30

  for _ = 1, max_iter do
    local p_close, weights, active_indices
    if #s == 1 then
      p_close = s[1]
      weights = { 1.0 }
      active_indices = { 1 }
    elseif #s == 2 then
      p_close, weights, active_indices = closest_point_segment(s[1], s[2])
    elseif #s == 3 then
      p_close, weights, active_indices = closest_point_triangle(s[1], s[2], s[3])
    end

    local dist = vec3.norm(p_close)

    -- Prune simplex to active features
    local new_s, new_sa, new_sb = {}, {}, {}
    local pt_a_close = vec3.zero()
    local pt_b_close = vec3.zero()
    for _, idx in ipairs(active_indices) do
      table.insert(new_s, s[idx])
      table.insert(new_sa, sa[idx])
      table.insert(new_sb, sb[idx])
    end
    for i = 1, #weights do
      pt_a_close = pt_a_close + sa[i] * weights[i]
      pt_b_close = pt_b_close + sb[i] * weights[i]
    end
    s, sa, sb = new_s, new_sa, new_sb

    if dist < tol then
      return 0.0, pt_a_close, pt_b_close, vec3.zero()
    end

    local normal = p_close * (1.0 / dist)
    local w, wa, wb = support_minkowski(shape_a, shape_b, -normal)
    local gap = dist - vec3.dot(w, normal)

    if gap < tol or math.abs(last_dist - dist) < 1e-9 then
      return dist, pt_a_close, pt_b_close, normal
    end

    last_dist = dist
    table.insert(s, w)
    table.insert(sa, wa)
    table.insert(sb, wb)

    if #s > 3 then
      table.remove(s, 1)
      table.remove(sa, 1)
      table.remove(sb, 1)
    end
  end

  local dist = vec3.norm(s[#s])
  return dist, sa[#sa], sb[#sb], (dist > 1e-12 and s[#s] * (1.0 / dist) or vec3.unit_z())
end

return gjk
