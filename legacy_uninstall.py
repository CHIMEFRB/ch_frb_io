#!/usr/bin/env python

import build_helpers

# If called recursively in superbuild, a global persistent LegacyUninstaller will be returned.
u = build_helpers.get_global_legacy_uninstaller()

u.uninstall_headers('assembled_chunk_msgpack.hpp')
u.uninstall_headers('ch_frb_io.hpp')
u.uninstall_headers('ch_frb_io_internals.hpp')
u.uninstall_headers('chlog.hpp')
u.uninstall_headers('L0_L1_packet.hpp')
u.uninstall_libraries('libch_frb_io*')
u.uninstall_executables('ch-show-intensity-file')
u.uninstall_executables('ch-plot-intensity-file')

# If called recursively in superbuild, run() will not be called here.
if __name__ == '__main__':
    lu.run()
