local csv_logger = {}
local logger_mt = { __index = csv_logger }

function csv_logger.new(headers)
  local l = {
    headers = headers or {},
    rows = {}
  }
  return setmetatable(l, logger_mt)
end

function csv_logger:log(row_data)
  local row = {}
  if #self.headers > 0 then
    for i, h in ipairs(self.headers) do
      local val = row_data[h] or row_data[i] or 0.0
      table.insert(row, val)
    end
  else
    for _, v in ipairs(row_data) do
      table.insert(row, v)
    end
  end
  table.insert(self.rows, row)
end

function csv_logger:to_csv_string()
  local lines = {}
  if #self.headers > 0 then
    table.insert(lines, table.concat(self.headers, ","))
  end

  for _, row in ipairs(self.rows) do
    local str_vals = {}
    for _, val in ipairs(row) do
      if type(val) == "number" then
        table.insert(str_vals, string.format("%.8g", val))
      else
        table.insert(str_vals, tostring(val))
      end
    end
    table.insert(lines, table.concat(str_vals, ","))
  end

  return table.concat(lines, "\n")
end

function csv_logger:write_file(filepath)
  local content = self:to_csv_string()
  local f, err = io.open(filepath, "w")
  if not f then
    error(string.format("Failed to open CSV file '%s': %s", tostring(filepath), tostring(err)))
  end
  f:write(content)
  f:close()
  return true
end

return csv_logger
