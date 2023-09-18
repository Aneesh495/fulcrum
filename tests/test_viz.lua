local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local mat3 = require("fulcrum.math.mat3")
local terminal = require("fulcrum.viz.terminal")
local rasterizer3d = require("fulcrum.viz.rasterizer3d")

runner.suite("Terminal Telemetry HUD and 3D Software Rasterizer")

runner.test("Terminal escape sequences and UI widgets", function()
  local move = terminal.move_to(5, 12)
  runner.assert_equal(move, "\27[5;12H", "ANSI cursor move sequence")

  local rgb = terminal.rgb_fg(255, 128, 0)
  runner.assert_equal(rgb, "\27[38;2;255;128;0m", "24-bit RGB truecolor sequence")

  local bar = terminal.progress_bar(0.5, 10)
  runner.assert_equal(bar, "█████░░░░░", "Progress bar 50% filled")

  local box = terminal.draw_box("TELEMETRY", { "Sim time: 1.25s", "Energy drift: 1e-7" }, 30)
  runner.assert_true(string.find(box, "TELEMETRY") ~= nil, "Box title rendered")
  runner.assert_true(string.find(box, "Sim time: 1.25s") ~= nil, "Box content line rendered")
end)

runner.test("3D Pinhole projection and camera clipping", function()
  local rast = rasterizer3d.new(80, 40)
  rast:set_camera(vec3.new(0, -5, 0), vec3.new(0, 0, 0), vec3.unit_z())

  -- Point at origin (straight ahead)
  local sx, sy, iz = rast:project(vec3.zero())
  runner.assert_true(sx ~= nil, "Point in front of camera projects")
  runner.assert_near(sx, 40.5, 1.5, "Center point x near midpoint of 80")
  runner.assert_near(sy, 20.5, 1.5, "Center point y near midpoint of 40")
  runner.assert_true(iz > 0.0, "Depth is positive")

  -- Point behind camera
  local behind_x = rast:project(vec3.new(0, -10, 0))
  runner.assert_true(behind_x == nil, "Point behind camera is clipped")
end)

runner.test("3D wireframe box rendering and depth buffer occlusion", function()
  local rast = rasterizer3d.new(60, 24)
  rast:set_camera(vec3.new(3, -4, 2), vec3.zero(), vec3.unit_z())

  rast:draw_grid(1.0, 0.5)
  rast:draw_axes(vec3.zero(), mat3.identity(), 0.5)
  rast:draw_box_wireframe(vec3.zero(), vec3.new(0.3, 0.3, 0.3), mat3.identity())

  local frame = rast:render_to_string()
  runner.assert_true(#frame > 500, "Frame buffer rendered to non-empty string")
  runner.assert_true(string.find(frame, "#") ~= nil, "Box wireframe vertices rendered")

  -- Test depth buffer overwrite
  rast:clear()
  rast:set_pixel(10, 10, 0.5, "A") -- farther
  rast:set_pixel(10, 10, 1.0, "B") -- closer (larger 1/z)
  runner.assert_equal(rast.char_buf[10][10], "B", "Closer pixel overwrote farther pixel")

  rast:set_pixel(10, 10, 0.2, "C") -- farther again
  runner.assert_equal(rast.char_buf[10][10], "B", "Farther pixel was rejected by z-buffer")
end)
