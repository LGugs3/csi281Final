cmake_minimum_required(VERSION 3.29)
project(csi281Final)

set(CMAKE_CXX_STANDARD 20)

#SET(FFMPEG_LIBRARY_PATH "C:/ffmpeg"/lib)
#INCLUDE_DIRECTORIES("C:/ffmpeg"/include)
#LINK_DIRECTORIES(${FFMPEG_LIBRARY_PATH})
#TARGET_LINK_LIBRARIES(${PROJECT_NAME}
#        ${FFMPEG_LIBRARY_PATH}/libavcodec.a
#        ${FFMPEG_LIBRARY_PATH}/libavutil.a
#        ${FFMPEG_LIBRARY_PATH}/libavdevice.a
#        ${FFMPEG_LIBRARY_PATH}/libswscale.a
#        ${FFMPEG_LIBRARY_PATH}/libswresample.a
#        pthread
#        z
#        m
#        dl
#)


add_executable(csi281Final
    main.cpp
        WAVFile.h
        WAVFile.cpp)