local terminal = {}

local ESC = "\27["

terminal.CLEAR_SCREEN = ESC .. "2J"
terminal.CLEAR_LINE   = ESC .. "2K"
terminal.CURSOR_HOME  = ESC .. "H"
terminal.CURSOR_HIDE  = ESC .. "?25l"
terminal.CURSOR_SHOW  = ESC .. "?25h"
terminal.RESET        = ESC .. "0m"
terminal.BOLD         = ESC .. "1m"
terminal.DIM          = ESC .. "2m"

function terminal.move_to(row, col)
  return string.format("%s%d;%dH", ESC, math.floor(row), math.floor(col))
end

function terminal.rgb_fg(r, g, b)
  return string.format("%s38;2;%d;%d;%dm", ESC, math.floor(r), math.floor(g), math.floor(b))
end

function terminal.rgb_bg(r, g, b)
  return string.format("%s48;2;%d;%d;%dm", ESC, math.floor(r), math.floor(g), math.floor(b))
end

function terminal.progress_bar(fraction, width)
  width = width or 20
  local filled = math.floor(math.max(0.0, math.min(1.0, fraction)) * width)
  local empty = width - filled
  return string.rep("█", filled) .. string.rep("░", empty)
end

function terminal.draw_box(title, lines, width)
  width = width or 50
  local out = {}
  local header = "┌─ " .. title .. " " .. string.rep("─", math.max(0, width - #title - 5)) .. "┐"
  table.insert(out, header)

  for _, line in ipairs(lines) do
    local stripped = string.gsub(line, "\27%[[%d;]*m", "")
    local pad = width - #stripped - 2
    if pad < 0 then pad = 0 end
    table.insert(out, "│ " .. line .. string.rep(" ", pad) .. " │")
  end

  local footer = "└" .. string.rep("─", width) .. "┘"
  table.insert(out, footer)
  return table.concat(out, "\n")
end

return terminal
