project(cudapoa_test)

set(sources
    main.cpp)

cuda_add_executable(cudapoa_test ${sources})

target_link_libraries(cudapoa_test cudapoa)
