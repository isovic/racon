.PHONY: all clean meson cmake debug dist modules

all: meson

clean:
	rm -rf build build-meson

meson:
	@echo "[Invoking Meson]"
	@mkdir -p build-meson && cd build-meson && meson --buildtype=release -Dc_args=-O3 && ninja

cmake:
	@echo "[Invoking CMake]"
	@mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

debug:
	@echo "[Invoking Meson]"
	@mkdir -p build-debug && cd build-debug && (meson --buildtype=debugoptimized -Db_sanitize=address) && ninja

dist: release
	cd build && ninja-dist

modules:
	git submodule update --init
