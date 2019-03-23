.PHONY: all clean release debug debug-gcc6 dist modules

all: release

clean:
	rm -rf build

release:
	@echo "[Invoking Meson]"
	@mkdir -p build && cd build && meson --buildtype=release -Dc_args=-O3 && ninja

debug:
	@echo "[Invoking Meson]"
	@mkdir -p build-debug && cd build-debug && (meson --buildtype=debugoptimized -Db_sanitize=address) && ninja

dist: release
	cd build && ninja-dist

modules:
	git submodule update --init
