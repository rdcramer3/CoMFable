cmake_minimum_required(VERSION 3.2)
project(dcomfa)

set (CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -m64 -W -Wall -stdlib=libc++ -dynamic -fPIC")
set(CMAKE_EXE_LINKER_FLAGS "-m64 -Wl,-rpath,@executable_path/. -Wl,-rpath,@executable_path/../../../")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")

set(SYSLIBS "-lpthread")
set(SYSZLIB "-lz")
set(SYSRASTERLIBS "-framework ApplicationServices")

set(SOURCE_FILES
    comfa.cpp)

set (OE_DIR "~/openeye/toolkits")

set (OE_LIBS
    oegrid
    oefizzchem
    oedepict
    oebio
    oechem
    oesystem
    oeplatform
)

include_directories ( ${OE_DIR}/include )
link_directories ( ${OE_DIR}/lib )

set (LIBS "${OE_LIBS};${SYSLIBS};${SYSRASTERLIBS};${SYSZLIB}")

set (SRCS
    "dcomfa ${SOURCE_FILES}"
)

include (${OE_DIR}/examples/ExamplesTemplate.cmake)

#add_executable(dcomfa ${SOURCE_FILES})
#target_link_libraries(dcomfa ${LIBS})
