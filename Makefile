CC = clang
DEBUG_FLAGS = -Wextra -Wall -fsanitize=address -g
RELEASE_FLAGS = -O3

THIRDPARTY_DIR = $(CURDIR)/thirdparty
RAYLIB_DIR = $(THIRDPARTY_DIR)/raylib-5.5_macos
GLFW_DIR = $(THIRDPARTY_DIR)/glfw-3.4_macos
RAYLIB_LIB_DIR = $(RAYLIB_DIR)/lib
GLFW_LIB_DIR = $(GLFW_DIR)/lib-arm64

LINK = -L$(RAYLIB_LIB_DIR) \
			 -L$(GLFW_LIB_DIR) \
			 -Wl,-rpath,$(RAYLIB_LIB_DIR)
LIBS = -lraylib -lglfw3
FRAMEWORKS = \
						 -framework Cocoa \
						 -framework OpenGL \
						 -framework IOKit \
						 -framework CoreVideo \
						 -framework QuartzCore


.PHONY = all main run clean debug sanitize check-deps

all: main run

main: check-deps
	mkdir -p bin
	$(CC) src/main.c -o bin/main $(LINK) $(LIBS) $(FRAMEWORKS) $(RELEASE_FLAGS)

debug: sanitize run

sanitize: check-deps
	mkdir -p bin
	$(CC) src/main.c -o bin/main $(LINK) $(LIBS) $(FRAMEWORKS) $(DEBUG_FLAGS)

run:
	./bin/main

clean:
	rm bin/*

check-deps:
	@test -d "$(RAYLIB_LIB_DIR)" || (echo "Missing raylib in $(RAYLIB_DIR). Run ./setup_macos.sh first." && exit 1)
	@test -d "$(GLFW_LIB_DIR)" || (echo "Missing glfw in $(GLFW_DIR). Run ./setup_macos.sh first." && exit 1)

