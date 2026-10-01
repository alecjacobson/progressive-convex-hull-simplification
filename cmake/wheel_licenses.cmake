# Include the notices for dependencies compiled into native Python wheels.
# scikit-build-core copies this metadata directory into the wheel's dist-info.
set(PCHS_LICENSE_DIR "${SKBUILD_METADATA_DIR}/licenses")

install(FILES "${libigl_SOURCE_DIR}/LICENSE.MPL2"
  DESTINATION "${PCHS_LICENSE_DIR}/libigl" COMPONENT pchs_python)
FetchContent_GetProperties(eigen SOURCE_DIR PCHS_EIGEN_SOURCE_DIR)
install(FILES
  "${PCHS_EIGEN_SOURCE_DIR}/COPYING.README"
  "${PCHS_EIGEN_SOURCE_DIR}/COPYING.MPL2"
  "${PCHS_EIGEN_SOURCE_DIR}/COPYING.BSD"
  "${PCHS_EIGEN_SOURCE_DIR}/COPYING.APACHE"
  "${PCHS_EIGEN_SOURCE_DIR}/COPYING.MINPACK"
  DESTINATION "${PCHS_LICENSE_DIR}/eigen" COMPONENT pchs_python)
install(FILES "${qhull_SOURCE_DIR}/COPYING.txt"
  DESTINATION "${PCHS_LICENSE_DIR}/qhull" COMPONENT pchs_python)
install(FILES "${nanobind_SOURCE_DIR}/LICENSE"
  DESTINATION "${PCHS_LICENSE_DIR}/nanobind" COMPONENT pchs_python)
install(FILES "${nanobind_SOURCE_DIR}/ext/robin_map/LICENSE"
  DESTINATION "${PCHS_LICENSE_DIR}/robin_map" COMPONENT pchs_python)

# These projects keep their notices in source files rather than LICENSE files.
install(FILES "${sdlp_SOURCE_DIR}/include/sdlp/sdlp.hpp"
  DESTINATION "${PCHS_LICENSE_DIR}/sdlp" COMPONENT pchs_python)
install(FILES "${predicates_SOURCE_DIR}/src/predicates.c"
  DESTINATION "${PCHS_LICENSE_DIR}/predicates" COMPONENT pchs_python)
