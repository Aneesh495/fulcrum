local models = {
  builder = require("fulcrum.models.builder"),
  urdf_parser = require("fulcrum.models.urdf_parser"),
  cartpole = require("fulcrum.models.cartpole"),
  pendulums = require("fulcrum.models.pendulums"),
  dzhanibekov = require("fulcrum.models.dzhanibekov"),
  newtons_cradle = require("fulcrum.models.newtons_cradle"),
  manipulator_7dof = require("fulcrum.models.manipulator_7dof"),
  biped = require("fulcrum.models.biped"),
  quadruped = require("fulcrum.models.quadruped")
}

return models
