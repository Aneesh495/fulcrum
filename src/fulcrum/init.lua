-- src/fulcrum/init.lua
-- Top-level module aggregator for the Fulcrum Geometric Mechanics and Multibody Dynamics Engine.

local fulcrum = {
  _VERSION     = "1.0.0",
  _DESCRIPTION = "Geometric mechanics, Lie group variational integrators, and articulated multibody dynamics",
  _AUTHOR      = "Fulcrum Contributors",
  _LICENSE     = "MIT"
}

-- Mathematical Foundations
fulcrum.math = {
  vec3       = require("fulcrum.math.vec3"),
  quat       = require("fulcrum.math.quat"),
  mat3       = require("fulcrum.math.mat3"),
  matrix     = require("fulcrum.math.matrix"),
  decomp     = require("fulcrum.math.decomp"),
  spatial    = require("fulcrum.math.spatial"),
  so3        = require("fulcrum.math.so3"),
  se3        = require("fulcrum.math.se3"),
  autodiff   = require("fulcrum.math.autodiff"),
  polynomial = require("fulcrum.math.polynomial")
}

-- Articulated Multibody Dynamics
fulcrum.dynamics = {
  model             = require("fulcrum.dynamics.model"),
  body              = require("fulcrum.dynamics.body"),
  joint             = require("fulcrum.dynamics.joint"),
  kinematics        = require("fulcrum.dynamics.kinematics"),
  jacobians         = require("fulcrum.dynamics.jacobians"),
  rnea              = require("fulcrum.dynamics.rnea"),
  crba              = require("fulcrum.dynamics.crba"),
  aba               = require("fulcrum.dynamics.aba"),
  operational_space = require("fulcrum.dynamics.operational_space")
}

-- Geometric & Symplectic Integrators
fulcrum.integrators = {
  symplectic     = require("fulcrum.integrators.symplectic"),
  forest_ruth    = require("fulcrum.integrators.forest_ruth"),
  gauss_legendre = require("fulcrum.integrators.gauss_legendre"),
  rkf45          = require("fulcrum.integrators.rkf45"),
  lgvi           = require("fulcrum.integrators.lie_group_integrator")
}

-- Non-Smooth Contact & Complementarity
fulcrum.contact = {
  friction    = require("fulcrum.contact.friction"),
  bilateral   = require("fulcrum.contact.bilateral"),
  lemke       = require("fulcrum.contact.lemke"),
  pgs         = require("fulcrum.contact.pgs"),
  moreau_jean = require("fulcrum.contact.moreau_jean")
}

-- Computational Geometry & Collision Detection
fulcrum.geometry = {
  shapes = require("fulcrum.geometry.shapes"),
  gjk    = require("fulcrum.geometry.gjk"),
  epa    = require("fulcrum.geometry.epa"),
  bvh    = require("fulcrum.geometry.bvh"),
  ccd    = require("fulcrum.geometry.ccd")
}

-- Optimal Control & Trajectory Optimization
fulcrum.opt = {
  ilqr        = require("fulcrum.opt.ilqr"),
  collocation = require("fulcrum.opt.collocation")
}

-- Articulated Robot Models & Parsers
fulcrum.models = {
  builder          = require("fulcrum.models.builder"),
  urdf_parser      = require("fulcrum.models.urdf_parser"),
  cartpole         = require("fulcrum.models.cartpole"),
  pendulums        = require("fulcrum.models.pendulums"),
  dzhanibekov      = require("fulcrum.models.dzhanibekov"),
  newtons_cradle   = require("fulcrum.models.newtons_cradle"),
  manipulator_7dof = require("fulcrum.models.manipulator_7dof"),
  biped            = require("fulcrum.models.biped"),
  quadruped        = require("fulcrum.models.quadruped")
}

-- Visualization & Telemetry
fulcrum.viz = {
  terminal     = require("fulcrum.viz.terminal"),
  rasterizer3d = require("fulcrum.viz.rasterizer3d"),
  obj_exporter = require("fulcrum.viz.obj_exporter"),
  csv_logger   = require("fulcrum.viz.csv_logger"),
  json_exporter= require("fulcrum.viz.json_exporter")
}

-- Native C Acceleration & LuaJIT FFI
fulcrum.ffi = {
  c_accelerator = require("fulcrum.ffi.c_accelerator")
}

return fulcrum
