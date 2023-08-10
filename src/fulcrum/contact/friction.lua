local vec3 = require("fulcrum.math.vec3")

local friction = {}

function friction.project_coulomb_cone(lambda_n, lambda_t1, lambda_t2, mu)
  -- Projects 3D contact impulse into circular Coulomb cone:
  -- K = { (fn, ft1, ft2) : fn >= 0, sqrt(ft1^2 + ft2^2) <= mu * fn }
  mu = mu or 0.5
  local fn = math.max(0.0, lambda_n)
  local ft_norm = math.sqrt(lambda_t1 * lambda_t1 + lambda_t2 * lambda_t2)

  local max_ft = mu * fn
  if ft_norm <= max_ft then
    return fn, lambda_t1, lambda_t2
  else
    if ft_norm > 1e-15 then
      local scale = max_ft / ft_norm
      return fn, lambda_t1 * scale, lambda_t2 * scale
    else
      return fn, 0.0, 0.0
    end
  end
end

function friction.linearized_cone_facets(normal, num_facets)
  -- Constructs tangential facet unit vectors d_k for polyhedral cone approximation
  num_facets = num_facets or 4
  local u_norm = vec3.unit(normal)

  -- Find tangent basis (t1, t2)
  local t1
  if math.abs(u_norm.x) < 0.9 then
    t1 = vec3.unit(vec3.cross(u_norm, vec3.unit_x()))
  else
    t1 = vec3.unit(vec3.cross(u_norm, vec3.unit_y()))
  end
  local t2 = vec3.cross(u_norm, t1)

  local facets = {}
  local d_theta = (2.0 * math.pi) / num_facets
  for k = 1, num_facets do
    local angle = (k - 1) * d_theta
    local c = math.cos(angle)
    local s = math.sin(angle)
    local d_k = (t1 * c) + (t2 * s)
    table.insert(facets, d_k)
  end

  return facets
end

return friction
