# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)
# 项目信息
project (ENU_K_BUNDLE)
# 查找当前目录下的所有源文件
# 并将名称保存到 DIR_SRCS 变量
aux_source_directory(. DIR_SRCS)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
# 指定生成目标
add_executable(main ${DIR_SRCS})
