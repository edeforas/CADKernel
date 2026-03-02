include_directories(../src)

add_executable(sample_nurbs_boolean sample_nurbs_boolean.cpp  )
target_link_libraries(sample_nurbs_boolean libCADTest)

add_executable(sample_nurbs_spheres_difference sample_nurbs_spheres_difference.cpp)
target_link_libraries(sample_nurbs_spheres_difference libCADTest)

add_executable(sample_bump_cube sample_bump_cube.cpp  )
target_link_libraries(sample_bump_cube libCADTest)

add_executable(sample_nurbs_mandelbrot sample_nurbs_mandelbrot.cpp )
target_link_libraries(sample_nurbs_mandelbrot libCADTest)

add_executable(sample_fractal_tetrahedron sample_fractal_tetrahedron.cpp )
target_link_libraries(sample_fractal_tetrahedron libCADTest)

add_executable(sample_mesh_mandelbrot sample_mesh_mandelbrot.cpp )
target_link_libraries(sample_mesh_mandelbrot libCADTest)

add_executable(sample_mesh_factory sample_mesh_factory.cpp )
target_link_libraries(sample_mesh_factory libCADTest)

add_executable(sample_nurbs_factory sample_nurbs_factory.cpp )
target_link_libraries(sample_nurbs_factory libCADTest)