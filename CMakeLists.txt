cmake_minimum_required(VERSION 3.16.0)
project (ctpmoran)

list(APPEND flags "-march=native")

set(SOURCE_FILES main.cpp
		 moran.cpp
                 state.cpp
                 population.cpp
                 strategy.cpp
                 ctpgame.cpp
                 ctpdata.cpp
                 rangen.cpp)

add_executable(ctpmoran ${SOURCE_FILES})

target_compile_options(ctpmoran
        PRIVATE
                ${flags}
)

find_package(GSL REQUIRED)    
target_link_libraries(ctpmoran GSL::gsl GSL::gslcblas)
