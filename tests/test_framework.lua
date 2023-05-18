local M = {}

local total_tests = 0
local passed_tests = 0
local failed_tests = 0
local suites = {}
local current_suite = nil

local ANSI_RESET  = "\27[0m"
local ANSI_RED    = "\27[31m"
local ANSI_GREEN  = "\27[32m"
local ANSI_YELLOW = "\27[33m"
local ANSI_CYAN   = "\27[36m"
local ANSI_BOLD   = "\27[1m"

function M.suite(name)
  current_suite = {
    name = name,
    tests = {}
  }
  table.insert(suites, current_suite)
  return current_suite
end

function M.test(name, fn)
  if not current_suite then
    M.suite("Default Suite")
  end
  table.insert(current_suite.tests, {
    name = name,
    fn = fn
  })
end

function M.assert_true(cond, msg)
  if not cond then
    error(msg or "Assertion failed: expected true condition", 2)
  end
end

function M.assert_false(cond, msg)
  if cond then
    error(msg or "Assertion failed: expected false condition", 2)
  end
end

function M.assert_equal(actual, expected, msg)
  if actual ~= expected then
    error(string.format("%s: expected %s, got %s",
      msg or "Equality assertion failed",
      tostring(expected),
      tostring(actual)), 2)
  end
end

function M.assert_near(actual, expected, tol, msg)
  tol = tol or 1e-7
  local diff = math.abs(actual - expected)
  if diff > tol then
    error(string.format("%s: expected %g, got %g (diff %g > tol %g)",
      msg or "Nearness assertion failed",
      expected, actual, diff, tol), 2)
  end
end

function M.assert_vector_near(actual, expected, tol, msg)
  tol = tol or 1e-7
  local n = #expected
  if #actual ~= n then
    error(string.format("%s: dimension mismatch (%d vs %d)",
      msg or "Vector dimension error", #actual, n), 2)
  end
  for i = 1, n do
    local diff = math.abs(actual[i] - expected[i])
    if diff > tol then
      error(string.format("%s at index %d: expected %g, got %g (diff %g > tol %g)",
        msg or "Vector nearness assertion failed",
        i, expected[i], actual[i], diff, tol), 2)
    end
  end
end

function M.run()
  print(string.format("%s%s=== Running Fulcrum Test Suite ===%s", ANSI_BOLD, ANSI_CYAN, ANSI_RESET))
  local start_time = os.clock()

  for _, s in ipairs(suites) do
    print(string.format("\n%s[%s]%s", ANSI_BOLD, s.name, ANSI_RESET))
    for _, t in ipairs(s.tests) do
      total_tests = total_tests + 1
      local ok, err = pcall(t.fn)
      if ok then
        passed_tests = passed_tests + 1
        print(string.format("  %s✓%s %s", ANSI_GREEN, ANSI_RESET, t.name))
      else
        failed_tests = failed_tests + 1
        print(string.format("  %s✗%s %s", ANSI_RED, ANSI_RESET, t.name))
        print(string.format("    %s%s%s", ANSI_RED, tostring(err), ANSI_RESET))
      end
    end
  end

  local elapsed = os.clock() - start_time
  print(string.format("\n%s----------------------------------------%s", ANSI_CYAN, ANSI_RESET))
  print(string.format("Total: %d | %sPassed: %d%s | %sFailed: %d%s | Time: %.4fs",
    total_tests,
    ANSI_GREEN, passed_tests, ANSI_RESET,
    failed_tests > 0 and ANSI_RED or ANSI_GREEN, failed_tests, ANSI_RESET,
    elapsed))

  if failed_tests > 0 then
    print(string.format("%sFAILED (%d tests failed)%s", ANSI_RED, failed_tests, ANSI_RESET))
    os.exit(1)
  else
    print(string.format("%sALL TESTS PASSED%s", ANSI_GREEN, ANSI_RESET))
  end
end

return M
