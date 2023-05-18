LUA ?= lua
LUAJIT ?= luajit
CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra -fPIC -std=c99

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SHLIB_EXT = dylib
	SHLIB_FLAGS = -dynamiclib -undefined dynamic_lookup
else
	SHLIB_EXT = so
	SHLIB_FLAGS = -shared
endif

C_SRC = src/c/fulcrum_c.c
C_OUT = src/c/libfulcrum.$(SHLIB_EXT)

.PHONY: all test test-jit bench clean c_lib help

all: c_lib test

c_lib: $(C_OUT)

$(C_OUT): $(C_SRC)
	@mkdir -p src/c
	$(CC) $(CFLAGS) $(SHLIB_FLAGS) -Isrc/c -o $@ $(C_SRC)

test:
	LUA_PATH="src/?.lua;src/?/init.lua;;" $(LUA) tests/run_tests.lua

test-jit: c_lib
	LUA_PATH="src/?.lua;src/?/init.lua;;" $(LUAJIT) tests/run_tests.lua

bench:
	LUA_PATH="src/?.lua;src/?/init.lua;;" $(LUA) bench/run_benchmarks.lua

clean:
	rm -f src/c/*.$(SHLIB_EXT) src/c/*.o
	rm -f *.csv *.obj *.jsonl

help:
	@echo "Fulcrum: Lie Group Variational Integrator and Multibody Dynamics Engine"
	@echo "Targets:"
	@echo "  make test      Run unit test suite using standard Lua"
	@echo "  make test-jit  Run unit test suite using LuaJIT and C accelerator"
	@echo "  make bench     Run performance benchmarks"
	@echo "  make c_lib     Compile optional C accelerator shared library"
	@echo "  make clean     Remove build artifacts and temporary files"
