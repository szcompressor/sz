#
# - Try to find Facebook zstd library
# This will define
# ZSTD_FOUND
# ZSTD_INCLUDE_DIRS
# ZSTD_LIBRARIES
#

find_path(
    ZSTD_INCLUDE_DIR
    NAMES "zstd.h"
    HINTS
        "/usr/local/facebook/include"
)

set (ZSTD_INCLUDE_DIRS ${ZSTD_INCLUDE_DIR})

find_library(
    ZSTD_LIBRARY
    NAMES zstd
    HINTS
        "/usr/local/facebook/lib"
)

set (ZSTD_LIBRARIES ${ZSTD_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    ZSTD ZSTD_INCLUDE_DIR ZSTD_LIBRARIES)

mark_as_advanced (ZSTD_INCLUDE_DIRS ZSTD_LIBRARIES)

