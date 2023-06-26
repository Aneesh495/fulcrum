local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local matrix = require("fulcrum.math.matrix")
local spatial = require("fulcrum.math.spatial")
local so3 = require("fulcrum.math.so3")
local se3 = require("fulcrum.math.se3")

local joint = {}

joint.TYPE_FIXED         = "fixed"
joint.TYPE_REVOLUTE      = "revolute"
joint.TYPE_PRISMATIC     = "prismatic"
joint.TYPE_HELICAL       = "helical"
joint.TYPE_SPHERICAL     = "spherical"
joint.TYPE_PLANAR        = "planar"
joint.TYPE_FLOATING_BASE = "floating_base"

function joint.revolute(axis, name)
  local u_axis = vec3.unit(axis or vec3.unit_z())
  return {
    type = joint.TYPE_REVOLUTE,
    name = name or "revolute_joint",
    dof = 1,
    axis = u_axis,
    limit_lower = -math.pi,
    limit_upper = math.pi,
    limit_velocity = 50.0,
    limit_effort = 1000.0,
    damping = 0.0,
    friction = 0.0,

    transform = function(q)
      local theta = (type(q) == "table") and q[1] or (q or 0.0)
      local R = mat3.from_axis_angle(u_axis, theta)
      return spatial.transform(R, vec3.zero())
    end,

    subspace = function(_)
      -- Spatial motion subspace matrix S in R^{6 x 1}
      -- S = [ axis^T, 0^T ]^T
      local S = matrix.zeros(6, 1)
      S[1][1] = u_axis.x
      S[2][1] = u_axis.y
      S[3][1] = u_axis.z
      return S
    end,

    bias_acceleration = function(_, _)
      -- c_J = S_dot * q_dot = 0 for single-axis revolute
      return spatial.zero()
    end
  }
end

function joint.prismatic(axis, name)
  local u_axis = vec3.unit(axis or vec3.unit_z())
  return {
    type = joint.TYPE_PRISMATIC,
    name = name or "prismatic_joint",
    dof = 1,
    axis = u_axis,
    limit_lower = -10.0,
    limit_upper = 10.0,
    limit_velocity = 20.0,
    limit_effort = 1000.0,
    damping = 0.0,
    friction = 0.0,

    transform = function(q)
      local d = (type(q) == "table") and q[1] or (q or 0.0)
      local p = u_axis * d
      return spatial.transform(mat3.identity(), p)
    end,

    subspace = function(_)
      -- S = [ 0^T, axis^T ]^T
      local S = matrix.zeros(6, 1)
      S[4][1] = u_axis.x
      S[5][1] = u_axis.y
      S[6][1] = u_axis.z
      return S
    end,

    bias_acceleration = function(_, _)
      return spatial.zero()
    end
  }
end

function joint.helical(axis, pitch, name)
  local u_axis = vec3.unit(axis or vec3.unit_z())
  local h = pitch or 0.01 -- translation per radian
  return {
    type = joint.TYPE_HELICAL,
    name = name or "helical_joint",
    dof = 1,
    axis = u_axis,
    pitch = h,
    limit_lower = -100.0,
    limit_upper = 100.0,
    limit_velocity = 50.0,
    limit_effort = 1000.0,

    transform = function(q)
      local theta = (type(q) == "table") and q[1] or (q or 0.0)
      local R = mat3.from_axis_angle(u_axis, theta)
      local p = u_axis * (h * theta)
      return spatial.transform(R, p)
    end,

    subspace = function(_)
      local S = matrix.zeros(6, 1)
      S[1][1] = u_axis.x
      S[2][1] = u_axis.y
      S[3][1] = u_axis.z
      S[4][1] = h * u_axis.x
      S[5][1] = h * u_axis.y
      S[6][1] = h * u_axis.z
      return S
    end,

    bias_acceleration = function(_, _)
      return spatial.zero()
    end
  }
end

function joint.spherical(name)
  return {
    type = joint.TYPE_SPHERICAL,
    name = name or "spherical_joint",
    dof = 3,
    limit_velocity = 50.0,
    limit_effort = 1000.0,

    transform = function(q)
      -- Parameterized by rotation vector q in R^3
      local phi = vec3.new(q[1] or 0.0, q[2] or 0.0, q[3] or 0.0)
      local R = so3.exp(phi)
      return spatial.transform(R, vec3.zero())
    end,

    subspace = function(q)
      local phi = vec3.new(q[1] or 0.0, q[2] or 0.0, q[3] or 0.0)
      local Jl = so3.left_jacobian(phi)
      local S = matrix.zeros(6, 3)
      for i = 1, 3 do
        for j = 1, 3 do
          S[i][j] = Jl[i][j]
        end
      end
      return S
    end,

    bias_acceleration = function(_, q_dot)
      -- Spatial velocity cross terms
      return spatial.zero()
    end
  }
end

function joint.floating_base(name)
  return {
    type = joint.TYPE_FLOATING_BASE,
    name = name or "floating_base_joint",
    dof = 6,

    transform = function(q)
      -- q = { px, py, pz, phi_x, phi_y, phi_z }
      local p = vec3.new(q[1] or 0.0, q[2] or 0.0, q[3] or 0.0)
      local phi = vec3.new(q[4] or 0.0, q[5] or 0.0, q[6] or 0.0)
      local R = so3.exp(phi)
      return spatial.transform(R, p)
    end,

    subspace = function(_)
      -- Free flyer: identity 6x6 motion subspace in body frame
      return matrix.eye(6)
    end,

    bias_acceleration = function(_, _)
      return spatial.zero()
    end
  }
end

function joint.fixed(name)
  return {
    type = joint.TYPE_FIXED,
    name = name or "fixed_joint",
    dof = 0,

    transform = function(_)
      return matrix.eye(6)
    end,

    subspace = function(_)
      return matrix.zeros(6, 0)
    end,

    bias_acceleration = function(_, _)
      return spatial.zero()
    end
  }
end

return joint
