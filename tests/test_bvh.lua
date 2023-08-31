local runner = require("tests.test_framework")
local vec3 = require("fulcrum.math.vec3")
local bvh = require("fulcrum.geometry.bvh")

runner.suite("Dynamic Bounding Volume Hierarchy (BVH)")

runner.test("AABB geometric properties and intersections", function()
  local b1 = bvh.aabb.new(vec3.new(0, 0, 0), vec3.new(2, 2, 2))
  local b2 = bvh.aabb.new(vec3.new(1, 1, 1), vec3.new(3, 3, 3))
  local b3 = bvh.aabb.new(vec3.new(4, 4, 4), vec3.new(5, 5, 5))

  runner.assert_near(bvh.aabb.surface_area(b1), 24.0, 1e-6, "Box surface area 2x(4+4+4)=24")
  runner.assert_true(bvh.aabb.intersects(b1, b2), "b1 and b2 intersect")
  runner.assert_false(bvh.aabb.intersects(b1, b3), "b1 and b3 do not intersect")

  local u = bvh.aabb.union(b1, b3)
  runner.assert_equal(u.min.x, 0.0, "Union min x")
  runner.assert_equal(u.max.x, 5.0, "Union max x")
  runner.assert_true(bvh.aabb.contains_point(u, vec3.new(2.5, 2.5, 2.5)), "Union contains midpoint")
end)

runner.test("Ray-AABB slab intersection", function()
  local box = bvh.aabb.new(vec3.new(1, -1, -1), vec3.new(3, 1, 1))
  local origin = vec3.new(0, 0, 0)
  local dir = vec3.new(1, 0, 0)

  local hit, tmin, tmax = bvh.aabb.ray_intersect(box, origin, dir)
  runner.assert_true(hit, "Ray hits box along +X")
  runner.assert_near(tmin, 1.0, 1e-6, "Entry distance tmin=1.0")
  runner.assert_near(tmax, 3.0, 1e-6, "Exit distance tmax=3.0")

  -- Ray pointing away
  local hit_away = bvh.aabb.ray_intersect(box, origin, vec3.new(-1, 0, 0))
  runner.assert_false(hit_away, "Ray pointing away misses box")
end)

runner.test("Tree insertion, query, and removal", function()
  local tree = bvh.new_tree()

  tree:insert("item_A", bvh.aabb.new(vec3.new(0, 0, 0), vec3.new(1, 1, 1)), { name = "A" })
  tree:insert("item_B", bvh.aabb.new(vec3.new(2, 0, 0), vec3.new(3, 1, 1)), { name = "B" })
  tree:insert("item_C", bvh.aabb.new(vec3.new(0.5, 0.5, 0), vec3.new(1.5, 1.5, 1)), { name = "C" })

  runner.assert_equal(tree.node_count, 3, "Node count is 3")

  -- Query overlapping box with A and C
  local query_box = bvh.aabb.new(vec3.new(0.2, 0.2, 0.2), vec3.new(0.8, 0.8, 0.8))
  local found = {}
  tree:query_aabb(query_box, function(id)
    found[id] = true
  end)

  runner.assert_true(found["item_A"], "Query box intersects item A")
  runner.assert_true(found["item_C"], "Query box intersects item C")
  runner.assert_false(found["item_B"], "Query box does not intersect item B")

  -- Broadphase pair query
  local pairs = {}
  tree:query_pairs(function(id1, id2)
    local k = id1 < id2 and (id1 .. ":" .. id2) or (id2 .. ":" .. id1)
    pairs[k] = true
  end)

  runner.assert_true(pairs["item_A:item_C"], "Pair A and C should collide")
  runner.assert_false(pairs["item_A:item_B"], "A and B should not collide")

  -- Raycast through tree
  local hits = tree:raycast(vec3.new(-1, 0.5, 0.5), vec3.new(1, 0, 0))
  runner.assert_true(#hits >= 2, "Ray hits at least 2 boxes along +X")
  runner.assert_equal(hits[1].user_id, "item_A", "Closest hit is item A")

  -- Removal
  local removed = tree:remove("item_A")
  runner.assert_true(removed, "Successfully removed item A")
  runner.assert_equal(tree.node_count, 2, "Node count decremented to 2")

  local hits_after = tree:raycast(vec3.new(-1, 0.5, 0.5), vec3.new(1, 0, 0))
  runner.assert_true(#hits_after >= 1, "Still hits remaining items")
  runner.assert_true(hits_after[1].user_id ~= "item_A", "Item A is no longer in hits")
end)
