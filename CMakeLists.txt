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


# # # 多项式测试程序
# # add_executable(test_ckks_poly 
# #     test/ckks_poly_test.c
# #     test/test_ckks_poly_main.c
# # )

# # target_include_directories(test_ckks_poly PUBLIC 
# #     ${CMAKE_CURRENT_SOURCE_DIR}/src
# # )
# # target_link_libraries(test_ckks_poly ckks)

# # 多项式测试程序 - 明确指定所有源文件
# add_executable(test_ckks_poly 
#     test/test_ckks_poly_main.c # 主函数
#     test/ckks_poly_test.c      # 测试实现
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

# 打印调试信息
message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")



# CKKS核心库源文件
file(GLOB CKKS_SOURCES "src/*.c")
message(STATUS "CKKS sources: ${CKKS_SOURCES}")

# Build encode library
add_library(ckks_encode STATIC src/ckks_encode.c)

# Build polynomial library (depends on encode)
add_library(ckks_polynomial STATIC src/ckks_polynomial.c)
target_link_libraries(ckks_polynomial ckks_encode m)
add_library(ckks_crypto STATIC src/ckks_crypto.c)
target_link_libraries(ckks_crypto ckks_polynomial ckks_encode)


# # Build main CKKS library (depends on polynomial and encode)
# add_library(ckks STATIC src/ckks_crypt.c src/ckks_keygen.c)
# target_link_libraries(ckks ckks_polynomial ckks_encode)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#本地openhitls构建的路径/usr/local
include_directories(SYSTEM /usr/local/include)

link_directories(/usr/local/lib)

add_executable(test_main test/test_main.c)

target_link_libraries(test_main ckks_crypto ckks_polynomial ckks_encode m)


target_link_libraries(test_main 
    hitls_crypto
    hitls_bsl
    boundscheck
)

#代码现在依赖 hitls_crypto 和 hitls_bsl，要将这些库链接到测试程序中。

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
add_executable(test_main test/test_main.c) # 占位符目标，防止空项目错误
target_link_libraries(test_main ckks_crypto ckks_polynomial ckks_encode  m)