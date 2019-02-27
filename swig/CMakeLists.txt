find_package(SWIG)
if(SWIG_FOUND)
  include(UseSWIG)


  option(BUILD_PYTHON_WRAPPER "build python wrapper" OFF)

  if(BUILD_PYTHON_WRAPPER)
    find_package(Python COMPONENTS Interpreter Development)
    set_property(SOURCE pysz.i PROPERTY CPLUSPLUS ON)

    swig_add_library(
      pysz
      LANGUAGE python
      SOURCES pysz.i pysz.cc
      )

    target_link_libraries(
      ${SWIG_MODULE_pysz_REAL_NAME} PUBLIC SZ Python::Python
      )
    target_include_directories(
      ${SWIG_MODULE_pysz_REAL_NAME} PRIVATE SZ Python::Python
      "${Python_SITELIB}/numpy/core/include/"
      ${CMAKE_CURRENT_SOURCE_DIR}
      )
  endif()
endif()
