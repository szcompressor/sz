find_package(SWIG)
if(SWIG_FOUND)
  include(UseSWIG)


  option(BUILD_PYTHON_WRAPPER "build python wrapper" OFF)

  if(BUILD_PYTHON_WRAPPER)
    find_package(Python COMPONENTS Interpreter Development NumPy)
    set_property(SOURCE pysz.i PROPERTY CPLUSPLUS ON)

    swig_add_library(
      pysz
      LANGUAGE python
      SOURCES pysz.i pysz.cc
      )

    target_link_libraries(
      pysz PUBLIC SZ Python::Python Python::NumPy
      )
    target_include_directories(
      pysz PRIVATE SZ Python::Python
      ${CMAKE_CURRENT_SOURCE_DIR}
      )

    get_property(swig_generated_module TARGET pysz PROPERTY SWIG_SUPPORT_FILES)
    install(TARGETS pysz DESTINATION ${Python_SITELIB})
    install(FILES ${swig_generated_module} DESTINATION ${Python_SITELIB})

  endif()
endif()
