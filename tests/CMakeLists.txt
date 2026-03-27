include_directories(../src)

add_executable(test_geometry test_geometry.cpp  )
target_link_libraries(test_geometry libCADTest)

add_executable(test_geometryutil test_geometryutil.cpp  )
target_link_libraries(test_geometryutil libCADTest)

add_executable(test_transform test_transform.cpp  )
target_link_libraries(test_transform libCADTest)

add_executable(test_smoothmeshfactory "test_smoothmeshfactory.cpp"  )
target_link_libraries(test_smoothmeshfactory libCADTest)

add_executable(test_mesh_fileformat "test_mesh_fileformat.cpp"  )
target_link_libraries(test_mesh_fileformat libCADTest)

add_executable(test_mesh_level_of_details test_mesh_level_of_details.cpp  )
target_link_libraries(test_mesh_level_of_details libCADTest)

add_executable(test_mesh_geosphere "test_mesh_geosphere.cpp"  )
target_link_libraries(test_mesh_geosphere libCADTest)

add_executable(test_meshoperations test_meshoperations.cpp  )
target_link_libraries(test_meshoperations libCADTest)

add_executable(test_meshboolean test_meshboolean.cpp  )
target_link_libraries(test_meshboolean libCADTest)

add_executable(test_measurements test_measurements.cpp  )
target_link_libraries(test_measurements libCADTest)

add_executable(test_renderer test_renderer.cpp )
target_link_libraries(test_renderer libCADTest)

add_executable(test_bmpfile test_bmpfile.cpp  )
target_link_libraries(test_bmpfile libCADTest)

add_executable(test_mesh_revolve test_mesh_revolve.cpp  )
target_link_libraries(test_mesh_revolve libCADTest)

add_executable(test_stepfile test_stepfile.cpp  )
target_link_libraries(test_stepfile libCADTest)

add_executable(test_quadbeziertriangle test_quadbeziertriangle.cpp  )
target_link_libraries(test_quadbeziertriangle libCADTest)

add_executable(test_nurbscurve test_nurbscurve.cpp  )
target_link_libraries(test_nurbscurve libCADTest)

add_executable(test_nurbsutil test_nurbsutil.cpp  )
target_link_libraries(test_nurbsutil libCADTest)

add_executable(test_nurbsutil_trimmedsurface test_nurbsutil_trimmedsurface.cpp  )
target_link_libraries(test_nurbsutil_trimmedsurface libCADTest)

add_executable(test_nurbsfactory test_nurbsfactory.cpp  )
target_link_libraries(test_nurbsfactory libCADTest)

add_executable(test_nurbssurface test_nurbssurface.cpp  )
target_link_libraries(test_nurbssurface libCADTest)

add_executable(test_nurbsextrude test_nurbsextrude.cpp  )
target_link_libraries(test_nurbsextrude libCADTest)

add_executable(test_nurbsruled test_nurbsruled.cpp  )
target_link_libraries(test_nurbsruled libCADTest)

add_executable(test_nurbsfillet test_nurbsfillet.cpp  )
target_link_libraries(test_nurbsfillet libCADTest)

add_executable(test_nurbsrevolve test_nurbsrevolve.cpp  )
target_link_libraries(test_nurbsrevolve libCADTest)

add_executable(test_nurbssweep test_nurbssweep.cpp  )
target_link_libraries(test_nurbssweep libCADTest)

add_executable(test_nurbsfitting test_nurbsfitting.cpp  )
target_link_libraries(test_nurbsfitting libCADTest)

add_executable(test_nurbsboolean test_nurbsboolean.cpp  )
target_link_libraries(test_nurbsboolean libCADTest)

add_executable(test_nurbs_to_bezier test_nurbs_to_bezier.cpp  )
target_link_libraries(test_nurbs_to_bezier libCADTest)

add_test(test_geometry test_geometry)
add_test(test_geometryutil test_geometryutil)
add_test(test_transform test_transform)
add_test(test_bmpfile test_bmpfile)
add_test(test_measurements test_measurements)
add_test(test_mesh_level_of_details test_mesh_level_of_details)
add_test(test_mesh_geosphere test_mesh_geosphere)
add_test(test_meshboolean test_meshboolean)
add_test(test_stepfile test_stepfile)
add_test(test_mesh_revolve test_mesh_revolve)

add_test(test_smoothmeshfactory test_smoothmeshfactory)
add_test(test_quadbeziertriangle test_quadbeziertriangle)

add_test(test_nurbscurve test_nurbscurve)
add_test(test_nurbsutil test_nurbsutil)
add_test(test_nurbsutil_trimmedsurface test_nurbsutil_trimmedsurface)
add_test(test_nurbfactory test_nurbsfactory)
add_test(test_nurbssurface test_nurbssurface)
add_test(test_nurbsextrude test_nurbsextrude)
add_test(test_nurbsruled test_nurbsruled)
add_test(test_nurbsfillet test_nurbsfillet)
add_test(test_nurbsrevolve test_nurbsrevolve)
add_test(test_nurbssweep test_nurbssweep)
add_test(test_nurbsfitting test_nurbsfitting)
add_test(test_nurbsboolean test_nurbsboolean)
add_test(test_nurbs_to_bezier test_nurbs_to_bezier)

add_test(test_renderer test_renderer)