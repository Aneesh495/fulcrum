local json_exporter = {}

local function serialize_value(val)
  local t = type(val)
  if t == "number" then
    if val ~= val then return "null" end -- NaN
    if val == math.huge or val == -math.huge then return "null" end
    return string.format("%.8g", val)
  elseif t == "string" then
    local escaped = string.gsub(val, '\\', '\\\\')
    escaped = string.gsub(escaped, '"', '\\"')
    escaped = string.gsub(escaped, '\n', '\\n')
    escaped = string.gsub(escaped, '\r', '\\r')
    escaped = string.gsub(escaped, '\t', '\\t')
    return '"' .. escaped .. '"'
  elseif t == "boolean" then
    return val and "true" or "false"
  elseif t == "table" then
    -- Check if table is sequential array
    local is_array = true
    local max_idx = 0
    local count = 0
    for k, _ in pairs(val) do
      count = count + 1
      if type(k) == "number" and math.floor(k) == k and k >= 1 then
        if k > max_idx then max_idx = k end
      else
        is_array = false
        break
      end
    end
    if is_array and max_idx == count then
      local parts = {}
      for i = 1, max_idx do
        table.insert(parts, serialize_value(val[i]))
      end
      return "[" .. table.concat(parts, ",") .. "]"
    else
      local parts = {}
      for k, v in pairs(val) do
        table.insert(parts, string.format('"%s":%s', tostring(k), serialize_value(v)))
      end
      return "{" .. table.concat(parts, ",") .. "}"
    end
  else
    return "null"
  end
end

function json_exporter.to_json(val)
  return serialize_value(val)
end

function json_exporter.export_trajectory(meta, frames)
  local doc = {
    metadata = meta or {},
    frames = frames or {}
  }
  return json_exporter.to_json(doc)
end

function json_exporter.write_file(filepath, data)
  local str = (type(data) == "string") and data or json_exporter.to_json(data)
  local f, err = io.open(filepath, "w")
  if not f then
    error(string.format("Failed to open JSON file '%s': %s", tostring(filepath), tostring(err)))
  end
  f:write(str)
  f:close()
  return true
end

return json_exporter
