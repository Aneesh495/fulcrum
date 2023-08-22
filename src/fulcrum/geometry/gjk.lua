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

return gjk
