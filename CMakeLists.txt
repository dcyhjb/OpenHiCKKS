# cmake_minimum_required(VERSION 3.13)  # Note: target_link_directories requires 3.13+
# project(CKKS VERSION 1.0.0 LANGUAGES C)

# set(CMAKE_C_STANDARD 11)
# set(CMAKE_C_STANDARD_REQUIRED ON)

# # CKKS source files
# set(CKKS_SOURCES
#     src/ckks_crypt.c
#     src/ckks_encode.c
#     src/ckks_keygen.c
#     src/ckks_poly.c
# )

# # Create CKKS library
# add_library(ckks STATIC ${CKKS_SOURCES})

# # Modern approach: target-specific include and link directories
# target_include_directories(ckks SYSTEM PUBLIC /usr/local/include)
# target_link_directories(ckks PUBLIC /usr/local/lib)

# # Link specific libraries (no need for full paths when using target_link_directories)
# target_link_libraries(ckks 
#     hitls_tls
#     hitls_crypto
#     # hitls_auth 
#     # hitls_pki
#     hitls_bsl
#     # boundscheck
# )


# # # Polynomial test program
# # add_executable(test_ckks_poly 
# #     test/ckks_poly_test.c
# #     test/test_ckks_poly_main.c
# # )

# # target_include_directories(test_ckks_poly PUBLIC 
# #     ${CMAKE_CURRENT_SOURCE_DIR}/src
# # )
# # target_link_libraries(test_ckks_poly ckks)

# # Polynomial test program - explicitly specify all source files
# add_executable(test_ckks_poly 
#     test/test_ckks_poly_main.c # Main function
#     test/ckks_poly_test.c      # Test implementation
# )
# target_include_directories(test_ckks_poly PUBLIC 
#     ${CMAKE_CURRENT_SOURCE_DIR}/src
#     ${CMAKE_CURRENT_SOURCE_DIR}/test
# )
# target_link_libraries(test_ckks_poly ckks)

# # # Test program
# # add_executable(ckks_poly_test, ckks_poly_test.c)

# # # add_executable(test_ckks test/test_ckks.c)
# # # target_link_libraries(test_ckks ckks)

cmake_minimum_required(VERSION 3.12)
project(CKKS VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Print debug information
message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")



# CKKS core library source files
file(GLOB CKKS_SOURCES "src/*.c")
message(STATUS "CKKS sources: ${CKKS_SOURCES}")

# Build encode library
add_library(ckks_encode STATIC src/ckks_encode.c)

# Build polynomial library (depends on encode)
add_library(ckks_polynomial STATIC src/ckks_polynomial.c)
target_link_libraries(ckks_polynomial ckks_encode m)
add_library(ckks_crypto STATIC src/ckks_crypto.c)
# target_include_directories(ckks_crypto 
#     PRIVATE 
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
#     /usr/local/include
#     /usr/local/include/hitls
# )
target_include_directories(ckks_crypto SYSTEM PUBLIC /usr/local/include /usr/local/include/hitls /usr/local/include/hitls/bsl /usr/local/include/hitls/crypto /usr/local/include/hitls/pki /usr/local/include/hitls/tls)
target_link_libraries(ckks_crypto ckks_polynomial ckks_encode)


# # Build main CKKS library (depends on polynomial and encode)
# add_library(ckks STATIC src/ckks_crypt.c src/ckks_keygen.c)
# target_link_libraries(ckks ckks_polynomial ckks_encode)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Local openhitls build path /usr/local
include_directories(SYSTEM /usr/local/include)
# include_directories(SYSTEM /usr/local/include/hitls)

link_directories(/usr/local/lib)

add_executable(test_main test/test_main.c)

target_link_libraries(test_main ckks_crypto ckks_polynomial ckks_encode m)

# Performance benchmark test program - evaluates performance metrics of CKKS operations
add_executable(test_benchmark test/test_benchmark.c)
target_link_libraries(test_benchmark ckks_crypto ckks_polynomial ckks_encode m)


target_link_libraries(test_main 
    hitls_crypto
    hitls_bsl
    boundscheck
)

target_link_libraries(test_benchmark 
    hitls_crypto
    hitls_bsl
    boundscheck
)

# Code now depends on hitls_crypto and hitls_bsl, these libraries need to be linked to the test programs.

# include_directories(/usr/local/include /usr/local/include/hitls /usr/local/include/hitls/bsl /usr/local/include/hitls/crypto /usr/local/include/hitls/pki /usr/local/include/hitls/tls)
# # Modern approach: target-specific include and link directories

# target_include_directories(ckks SYSTEM PUBLIC /usr/local/include /usr/local/include/hitls /usr/local/include/hitls/bsl /usr/local/include/hitls/crypto /usr/local/include/hitls/pki /usr/local/include/hitls/tls)
# link_directories(/usr/local/lib)
# link_libraries(hitls_auth
#   hitls_tls
#     hitls_crypto
#     hitls_pki
#     hitls_bsl
#     boundscheck
# )
# target_link_libraries(ckks PUBLIC     
#     hitls_tls
#     hitls_crypto
#     hitls_auth
#     hitls_pki
#     hitls_bsl
#     boundscheck
# )



# link_libraries(ckks)
# add_executable(test_main test/test_main.c) # Placeholder target to prevent empty project error
# target_link_libraries(test_main ckks_crypto ckks_polynomial ckks_encode  m)