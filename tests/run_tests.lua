local runner = require("tests.test_framework")

-- Core math suites
require("tests.test_vec3")
require("tests.test_quat")
require("tests.test_mat3")
require("tests.test_matrix")
require("tests.test_decomp")
require("tests.test_spatial")
require("tests.test_so3")
require("tests.test_se3")
require("tests.test_autodiff")
require("tests.test_polynomial")

-- Dynamics suites
require("tests.test_joint")
require("tests.test_kinematics")
require("tests.test_jacobians")

runner.suite("Core Framework Initialization")
runner.test("Version metadata check", function()
  local ver = require("fulcrum.version")
  runner.assert_equal(ver.major, 1, "Major version")
  runner.assert_equal(ver.minor, 0, "Minor version")
  runner.assert_equal(ver.patch, 0, "Patch version")
  runner.assert_equal(ver.string, "1.0.0", "Version string")
  runner.assert_equal(ver.codename, "Lagrange", "Version codename")
end)

runner.run()
