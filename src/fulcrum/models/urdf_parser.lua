local vec3 = require("fulcrum.math.vec3")
local builder = require("fulcrum.models.builder")

local urdf_parser = {}

local function parse_numbers(str, count)
  if not str then return nil end
  local nums = {}
  for num in string.gmatch(str, "%S+") do
    local val = tonumber(num)
    if val then
      table.insert(nums, val)
    end
  end
  if count and #nums ~= count then
    return nil
  end
  return nums
end

local function parse_vec3_str(str)
  local nums = parse_numbers(str, 3)
  if nums then
    return vec3.new(nums[1], nums[2], nums[3])
  end
  return vec3.zero()
end

local function extract_attr(tag_str, attr_name)
  local pat1 = attr_name .. '%s*=%s*"([^"]*)"'
  local val = string.match(tag_str, pat1)
  if val then return val end
  local pat2 = attr_name .. "%s*=%s*'([^']*)'"
  return string.match(tag_str, pat2)
end

function urdf_parser.parse_string(xml_text)
  -- Extract robot name
  local robot_name = extract_attr(string.match(xml_text, "<robot%s+[^>]*>") or "", "name") or "robot"

  local links = {}
  local joints = {}

  -- Extract all <link> blocks
  for link_block in string.gmatch(xml_text, "<link%s+([^>]*)>(.-)</link>") do
    local link_header, link_body = link_block, ""
  end

  -- Comprehensive link extraction supporting both self-closing and paired tags
  for full_link in string.gmatch(xml_text, "<link(.-)</link>") do
    local header = string.match(full_link, "^(.-)>") or ""
    local name = extract_attr(header, "name")

    if name then
      local mass = 0.0
      local com = vec3.zero()
      local inertia = { ixx = 0.1, iyy = 0.1, izz = 0.1, ixy = 0, ixz = 0, iyz = 0 }

      local inertial_block = string.match(full_link, "<inertial>(.-)</inertial>")
      if inertial_block then
        local m_str = extract_attr(inertial_block, "value")
        if m_str then mass = tonumber(m_str) or 0.0 end

        local orig_tag = string.match(inertial_block, "<origin%s+[^>]*/>")
        if orig_tag then
          com = parse_vec3_str(extract_attr(orig_tag, "xyz"))
        end

        local in_tag = string.match(inertial_block, "<inertia%s+[^>]*/>")
        if in_tag then
          inertia.ixx = tonumber(extract_attr(in_tag, "ixx")) or 0.1
          inertia.iyy = tonumber(extract_attr(in_tag, "iyy")) or 0.1
          inertia.izz = tonumber(extract_attr(in_tag, "izz")) or 0.1
          inertia.ixy = tonumber(extract_attr(in_tag, "ixy")) or 0.0
          inertia.ixz = tonumber(extract_attr(in_tag, "ixz")) or 0.0
          inertia.iyz = tonumber(extract_attr(in_tag, "iyz")) or 0.0
        end
      end

      links[name] = {
        name = name,
        mass = mass,
        com = com,
        inertia = inertia
      }
    end
  end

  -- Extract all <joint> blocks
  for full_joint in string.gmatch(xml_text, "<joint(.-)</joint>") do
    local header = string.match(full_joint, "^(.-)>") or ""
    local name = extract_attr(header, "name")
    local j_type = extract_attr(header, "type") or "revolute"

    if name then
      local parent_tag = string.match(full_joint, "<parent%s+[^>]*/>")
      local child_tag = string.match(full_joint, "<child%s+[^>]*/>")
      local origin_tag = string.match(full_joint, "<origin%s+[^>]*/>")
      local axis_tag = string.match(full_joint, "<axis%s+[^>]*/>")
      local limit_tag = string.match(full_joint, "<limit%s+[^>]*/>")

      local parent_link = parent_tag and extract_attr(parent_tag, "link") or nil
      local child_link = child_tag and extract_attr(child_tag, "link") or nil

      local origin_xyz = vec3.zero()
      local origin_rpy = vec3.zero()
      if origin_tag then
        origin_xyz = parse_vec3_str(extract_attr(origin_tag, "xyz"))
        origin_rpy = parse_vec3_str(extract_attr(origin_tag, "rpy"))
      end

      local axis = vec3.unit_z()
      if axis_tag then
        axis = parse_vec3_str(extract_attr(axis_tag, "xyz"))
      end

      local limits = nil
      if limit_tag then
        limits = {
          lower = tonumber(extract_attr(limit_tag, "lower")),
          upper = tonumber(extract_attr(limit_tag, "upper")),
          effort = tonumber(extract_attr(limit_tag, "effort")),
          velocity = tonumber(extract_attr(limit_tag, "velocity"))
        }
      end

      table.insert(joints, {
        name = name,
        type = j_type,
        parent = parent_link,
        child = child_link,
        xyz = origin_xyz,
        rpy = origin_rpy,
        axis = axis,
        limits = limits
      })
    end
  end

  -- Find root link: link that is never a child of any joint
  local is_child = {}
  for _, j in ipairs(joints) do
    if j.child then is_child[j.child] = true end
  end

  local root_name = nil
  for name, _ in pairs(links) do
    if not is_child[name] then
      root_name = name
      break
    end
  end

  if not root_name then
    -- Fallback: first joint parent or default base_link
    root_name = joints[1] and joints[1].parent or "base_link"
  end

  -- Construct model using builder
  local b = builder.new(robot_name)
  local root_info = links[root_name] or { name = root_name, mass = 0.0, com = vec3.zero() }
  b:add_base(root_name, {
    mass = root_info.mass,
    com = root_info.com,
    inertia = root_info.inertia
  })

  -- Map child link to its incoming joint
  local child_to_joint = {}
  for _, j in ipairs(joints) do
    if j.child then child_to_joint[j.child] = j end
  end

  -- Add links in tree hierarchy
  for name, lk in pairs(links) do
    if name ~= root_name then
      local j = child_to_joint[name]
      if j then
        b:add_link({
          name = name,
          parent = j.parent,
          joint_name = j.name,
          joint_type = j.type,
          axis = j.axis,
          origin = { xyz = j.xyz, rpy = j.rpy },
          mass = lk.mass,
          com = lk.com,
          inertia = lk.inertia,
          limits = j.limits
        })
      end
    end
  end

  local model, meta = b:build()
  meta.joints = joints
  meta.links = links
  return model, meta
end

function urdf_parser.parse_file(filepath)
  local f, err = io.open(filepath, "r")
  if not f then
    error(string.format("Failed to open URDF file '%s': %s", tostring(filepath), tostring(err)))
  end
  local content = f:read("*a")
  f:close()
  return urdf_parser.parse_string(content)
end

return urdf_parser
