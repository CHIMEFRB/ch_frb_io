#!/usr/bin/env python

import build_helpers

# If called recursively in superbuild, a global persistent HeavyHandedUninstaller will be returned.
u = build_helpers.get_global_heavy_handed_uninstaller()

u.uninstall_headers('assembled_chunk_msgpack.hpp')
u.uninstall_headers('ch_frb_io.hpp')
u.uninstall_headers('ch_frb_io_internals.hpp')
u.uninstall_headers('chlog.hpp')
u.uninstall_headers('L0_L1_packet.hpp')
u.uninstall_headers('bitshuffle*.h', comment='heavy-handed uninstaller assumes that bitshuffle installation was done by ch_frb_io')
u.uninstall_headers('iochain.h', comment='heavy-handed uninstaller assumes that bitshuffle installation was done by ch_frb_io')
u.uninstall_libraries('libch_frb_io*')
u.uninstall_executables('ch-show-intensity-file')
u.uninstall_executables('ch-plot-intensity-file')

# If called recursively in superbuild, run() will not be called here.
if __name__ == '__main__':
    u.run()
