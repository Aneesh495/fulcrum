local integrators = {
  symplectic = require("fulcrum.integrators.symplectic"),
  forest_ruth = require("fulcrum.integrators.forest_ruth"),
  lie_group_integrator = require("fulcrum.integrators.lie_group_integrator"),
  gauss_legendre = require("fulcrum.integrators.gauss_legendre"),
  rkf45 = require("fulcrum.integrators.rkf45")
}

return integrators
