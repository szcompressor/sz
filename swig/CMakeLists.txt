find_package(SWIG)
if(SWIG_FOUND)
  include(UseSWIG)

  option(BUILD_PYTHON_WRAPPER "build python wrapper" OFF)

  if(BUILD_PYTHON_WRAPPER)
    find_package(Python3 COMPONENTS Interpreter Development)
    set_property(SOURCE pysz.i PROPERTY CPLUSPLUS ON)

    swig_add_library(
      pysz
      LANGUAGE python
      SOURCES pysz.i
      )

    target_link_libraries(
      ${SWIG_MODULE_pysz_REAL_NAME} PUBLIC SZ Python3::Python
      )
    target_include_directories(
      ${SWIG_MODULE_pysz_REAL_NAME} PRIVATE SZ Python3::Python
      )
  endif()
endif()
