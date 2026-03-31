#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THIRDPARTY_DIR="$SCRIPT_DIR/thirdparty"

RAYLIB_VERSION="5.5"
GLFW_VERSION="3.4"

RAYLIB_ARCHIVE="raylib-${RAYLIB_VERSION}_macos.tar.gz"
RAYLIB_URL="https://github.com/raysan5/raylib/releases/download/${RAYLIB_VERSION}/${RAYLIB_ARCHIVE}"
RAYLIB_DIR="$THIRDPARTY_DIR/raylib-${RAYLIB_VERSION}_macos"

GLFW_ARCHIVE="glfw-${GLFW_VERSION}.bin.MACOS.zip"
GLFW_URL="https://github.com/glfw/glfw/releases/download/${GLFW_VERSION}/${GLFW_ARCHIVE}"
GLFW_EXTRACTED_DIR="$THIRDPARTY_DIR/glfw-${GLFW_VERSION}.bin.MACOS"
GLFW_TARGET_DIR="$THIRDPARTY_DIR/glfw-${GLFW_VERSION}_macos"

mkdir -p "$THIRDPARTY_DIR"

download_file() {
  local url="$1"
  local out="$2"
  if [[ -f "$out" ]]; then
    echo "Found $(basename "$out"), skipping download"
    return
  fi

  echo "Downloading $(basename "$out")..."
  curl -L --fail --progress-bar "$url" -o "$out"
}

extract_raylib() {
  local archive_path="$THIRDPARTY_DIR/$RAYLIB_ARCHIVE"
  if [[ -d "$RAYLIB_DIR" ]]; then
    echo "Found $(basename "$RAYLIB_DIR"), skipping extract"
    return
  fi

  echo "Extracting $RAYLIB_ARCHIVE..."
  tar -xzf "$archive_path" -C "$THIRDPARTY_DIR"
}

extract_glfw() {
  local archive_path="$THIRDPARTY_DIR/$GLFW_ARCHIVE"

  if [[ -d "$GLFW_TARGET_DIR" ]]; then
    echo "Found $(basename "$GLFW_TARGET_DIR"), skipping extract"
    return
  fi

  echo "Extracting $GLFW_ARCHIVE..."
  unzip -q -o "$archive_path" -d "$THIRDPARTY_DIR"

  if [[ -d "$GLFW_EXTRACTED_DIR" ]]; then
    mv "$GLFW_EXTRACTED_DIR" "$GLFW_TARGET_DIR"
  fi
}

verify_deps() {
  [[ -d "$RAYLIB_DIR/lib" ]] || { echo "raylib lib directory not found"; exit 1; }
  [[ -d "$GLFW_TARGET_DIR/lib-arm64" ]] || { echo "glfw arm64 lib directory not found"; exit 1; }
}

download_file "$RAYLIB_URL" "$THIRDPARTY_DIR/$RAYLIB_ARCHIVE"
extract_raylib

download_file "$GLFW_URL" "$THIRDPARTY_DIR/$GLFW_ARCHIVE"
extract_glfw

verify_deps
echo "Setup complete. Run: make"
