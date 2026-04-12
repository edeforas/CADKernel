include_directories(../src)

add_executable(sample_nurbs_boolean sample_nurbs_boolean.cpp  )
target_link_libraries(sample_nurbs_boolean libCADTest)

add_executable(sample_nurbs_spheres_difference sample_nurbs_spheres_difference.cpp)
target_link_libraries(sample_nurbs_spheres_difference libCADTest)

add_executable(sample_bump_cube sample_bump_cube.cpp  )
target_link_libraries(sample_bump_cube libCADTest)

add_executable(sample_nurbs_mandelbrot sample_nurbs_mandelbrot.cpp )
target_link_libraries(sample_nurbs_mandelbrot libCADTest)

add_executable(sample_nurbs_offseting sample_nurbs_offseting.cpp )
target_link_libraries(sample_nurbs_offseting libCADTest)

add_executable(sample_fractal_tetrahedron sample_fractal_tetrahedron.cpp )
target_link_libraries(sample_fractal_tetrahedron libCADTest)

add_executable(sample_mesh_mandelbrot sample_mesh_mandelbrot.cpp )
target_link_libraries(sample_mesh_mandelbrot libCADTest)

add_executable(sample_mesh_factory sample_mesh_factory.cpp )
target_link_libraries(sample_mesh_factory libCADTest)

add_executable(sample_nurbs_factory sample_nurbs_factory.cpp )
target_link_libraries(sample_nurbs_factory libCADTest)

add_executable(sample_nurbs_fitting sample_nurbs_fitting.cpp )
target_link_libraries(sample_nurbs_fitting libCADTest)

add_executable(sample_nurbs_surface_degree sample_nurbs_surface_degree.cpp )
target_link_libraries(sample_nurbs_surface_degree libCADTest)

add_executable(sample_nurbs_sweep sample_nurbs_sweep.cpp )
target_link_libraries(sample_nurbs_sweep libCADTest)

add_executable(sample_nurbs_fillet sample_nurbs_fillet.cpp )
target_link_libraries(sample_nurbs_fillet libCADTest)

add_executable(sample_nurbs_extrude sample_nurbs_extrude.cpp )
target_link_libraries(sample_nurbs_extrude libCADTest)

add_executable(sample_mobius_strip sample_mobius_strip.cpp )
target_link_libraries(sample_mobius_strip libCADTest)

add_executable(sample_nurbs_cube_cylinder_difference sample_nurbs_cube_cylinder_difference.cpp)
target_link_libraries(sample_nurbs_cube_cylinder_difference libCADTest)

add_executable(sample_nurbs_curve_continuity sample_nurbs_curve_continuity.cpp)
target_link_libraries(sample_nurbs_curve_continuity libCADTest)