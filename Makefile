.PHONY: all clean meson cmake debug dist modules

all: meson

clean:
	rm -rf build build-meson

meson: modules
	@echo "[Invoking Meson]"
	@mkdir -p build-meson && cd build-meson && meson --buildtype=release -Dc_args=-O3 && ninja

rebuild: modules
	@echo "[Running Ninja only]"
	@ninja -C build-meson

cmake: modules
	@echo "[Invoking CMake]"
	@mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

debug: modules
	@echo "[Invoking Meson]"
	@mkdir -p build-debug && cd build-debug && (meson --buildtype=debugoptimized -Db_sanitize=address) && ninja

dist: release
	cd build && ninja-dist

modules:
	@echo "[Fetching submodules]"
	@git submodule update --init
