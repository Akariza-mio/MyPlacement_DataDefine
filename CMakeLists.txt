cmake_minimum_required(VERSION 3.16)
project(myplacement LANGUAGES CXX)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# 当前阶段只有命令行布局，不需要链接尚未使用的 X11 绘图库。
add_executable(myplace
    myplace_main.cpp
    MyPlace/myplace.cpp
    MyPlace/initial_place_io.cpp
    MyPlace/placement_report.cpp
    PlaceData/placedata.cpp
)
target_include_directories(myplace PRIVATE MyPlace PlaceData Common)
target_link_libraries(myplace PRIVATE Eigen3::Eigen)
if(WIN32)
    # 读取原始Unicode命令行参数，避免文件名在进入程序前被替换。
    target_link_libraries(myplace PRIVATE shell32)
endif()
set_target_properties(myplace PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
if(MSVC)
    target_compile_options(myplace PRIVATE /W4 /utf-8)
else()
    target_compile_options(myplace PRIVATE -Wall -Wextra -pedantic)
endif()

# test_data 按团队约定被忽略；正式构建不引用仅存在于本机的测试文件。
