local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local obj_exporter = require("fulcrum.viz.obj_exporter")
local csv_logger = require("fulcrum.viz.csv_logger")
local json_exporter = require("fulcrum.viz.json_exporter")

runner.suite("Data and Mesh Exporters (OBJ, CSV, JSON)")

runner.test("Wavefront OBJ 3D box mesh generation", function()
  local verts, faces = obj_exporter.box_mesh(vec3.new(1.0, 1.0, 1.0))
  runner.assert_equal(#verts, 8, "Box has 8 vertices")
  runner.assert_equal(#faces, 12, "Box has 12 triangular faces")

  local obj_str = obj_exporter.to_obj_string({
    { name = "cube", vertices = verts, faces = faces }
  })

  runner.assert_true(string.find(obj_str, "o cube") ~= nil, "Contains object header")
  runner.assert_true(string.find(obj_str, "v 1.000000 1.000000 1.000000") ~= nil, "Contains vertex coordinate")
  runner.assert_true(string.find(obj_str, "f 1 2 3") ~= nil, "Contains face definition")
end)

runner.test("CSV trajectory telemetry logger", function()
  local logger = csv_logger.new({ "time", "q1", "q2", "energy" })
  logger:log({ time = 0.0, q1 = 0.1, q2 = -0.2, energy = 4.5 })
  logger:log({ time = 0.01, q1 = 0.105, q2 = -0.198, energy = 4.50001 })

  local csv_str = logger:to_csv_string()
  runner.assert_true(string.find(csv_str, "time,q1,q2,energy", 1, true) ~= nil, "Contains CSV header")
  runner.assert_true(string.find(csv_str, "0.01,0.105,-0.198", 1, true) ~= nil, "Contains logged row values")
end)

runner.test("JSON serializer and trajectory format", function()
  local data = {
    robot = "arm_3dof",
    dof = 3,
    active = true,
    q = { 0.1, 0.2, -0.3 }
  }

  local json_str = json_exporter.to_json(data)
  runner.assert_true(string.find(json_str, '"robot":"arm_3dof"') ~= nil, "String property serialized")
  runner.assert_true(string.find(json_str, '"dof":3') ~= nil, "Number property serialized")
  runner.assert_true(string.find(json_str, '"active":true') ~= nil, "Boolean property serialized")
  runner.assert_true(string.find(json_str, '[0.1,0.2,-0.3]') ~= nil, "Array property serialized")

  local traj_json = json_exporter.export_trajectory({ name = "test_sim" }, { { t = 0.0 }, { t = 0.1 } })
  runner.assert_true(string.find(traj_json, '"metadata"') ~= nil, "Trajectory metadata tag")
  runner.assert_true(string.find(traj_json, '"frames"') ~= nil, "Trajectory frames array")
end)
