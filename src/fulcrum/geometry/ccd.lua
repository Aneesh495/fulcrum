local vec3 = require("fulcrum.math.vec3")
local shapes = require("fulcrum.geometry.shapes")
local gjk = require("fulcrum.geometry.gjk")

local ccd = {}

function ccd.solve_linear(shape_a, v_a, shape_b, v_b, dt, tol, max_iter)
  dt = dt or 1.0
  tol = tol or 1e-4
  max_iter = max_iter or 50

  local v_rel = v_b - v_a
  local v_max = vec3.norm(v_rel)

  if v_max < 1e-12 then
    -- No relative motion: check static overlap
    local hit = gjk.intersect(shape_a, shape_b)
    if hit then
      return true, 0.0, vec3.unit_z(), vec3.zero()
    else
      return false, dt, nil, nil
    end
  end

  local t = 0.0
  local current_a = shape_a
  local current_b = shape_b

  for _ = 1, max_iter do
    -- Query distance between shapes at current time t
    local dist, pt_a, pt_b, normal = gjk.distance(current_a, current_b, 40, tol * 0.1)

    if dist <= tol then
      -- Impact detected at current time t
      local contact_pt = pt_a and (pt_a + pt_b) * 0.5 or vec3.zero()
      return true, t, normal, contact_pt
    end

    -- Conservative advancement step: advance by distance / approach_speed
    local approach_speed = v_max
    if normal then
      local n_dot_v = vec3.dot(v_rel, normal)
      if n_dot_v < 0 then
        approach_speed = math.min(v_max, -n_dot_v)
      end
    end

    if approach_speed < 1e-12 then
      -- Shapes are separating
      return false, dt, nil, nil
    end

    local step = 0.9 * (dist / approach_speed)
    step = math.max(step, tol * 0.1)
    t = t + step

    if t > dt then
      -- Traversed full time window with no impact
      return false, dt, nil, nil
    end

    -- Advance shapes along trajectories
    current_a = shapes.translate(shape_a, v_a * t)
    current_b = shapes.translate(shape_b, v_b * t)
  end

  -- Max iterations reached without confirmed impact
  return false, dt, nil, nil
end

function ccd.solve_general(trajectory_a_fn, trajectory_b_fn, v_bound, dt, tol, max_iter)
  dt = dt or 1.0
  tol = tol or 1e-4
  max_iter = max_iter or 50
  v_bound = v_bound or 1.0

  local t = 0.0

  for _ = 1, max_iter do
    local s_a = trajectory_a_fn(t)
    local s_b = trajectory_b_fn(t)

    local dist, pt_a, pt_b, normal = gjk.distance(s_a, s_b, 40, tol * 0.1)
    if dist <= tol then
      local pt = pt_a and (pt_a + pt_b) * 0.5 or vec3.zero()
      return true, t, normal, pt
    end

    local step = 0.85 * (dist / v_bound)
    step = math.max(step, tol * 0.1)
    t = t + step

    if t > dt then
      return false, dt, nil, nil
    end
  end

  return false, dt, nil, nil
end

return ccd
