set (SOURCES time0r.cpp)
set (TARGET time0r)

set (CXX_FLAGS $(CXX_FLAGS) -O3)

if (MSVC)
  set (SOURCES ${SOURCES} ${FREI0R_DEF})
endif (MSVC)

add_library (${TARGET}  MODULE ${SOURCES})

# No «lib» prefix (name.so instead of libname.so)
set_target_properties (${TARGET} PROPERTIES PREFIX "")

install (TARGETS ${TARGET} LIBRARY DESTINATION ${LIBDIR})
