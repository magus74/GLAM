#include <scalar_field.hpp>


namespace vic {


  void scalar_field::import_from(const std_triangle_mesh& in_mesh) {
    mesh_.clear();
    std::size_t max_id = 0;
    for( std::size_t  i = 0; i < in_mesh.triangle_count(); ++i) {
      std_triangle_mesh::tri_t t  = in_mesh.get_triangle(i);
      mesh_.insert_triangle(tmesh_t::triangle_t(t[0],t[1],t[2]));
      for( std::size_t v = 0; v < 3; ++v) {
	if ( t[v] > max_id ) max_id = t[v];
	tmesh_t::vertex_data_t* v_ptr =  mesh_.vertex_data(t[v]);
	v_ptr->position() = in_mesh.get_position(t[v]);
	std_triangle_mesh::color_t c = in_mesh.get_color(t[v]);
	v_ptr->color() = tmesh_t::vertex_data_t::color_t(c[0],c[1],c[2]);
	v_ptr->normal() = in_mesh.get_normal(t[v]);
      }
    }
    std::cerr<< "Max id  == " << max_id << std::endl;
    std::cerr<< "Imported mesh with " << mesh_.triangle_count() << " triangles " << std::endl;    

  }

  void scalar_field::export_to(std_triangle_mesh& out_mesh) {

    std::size_t nfaces = mesh_.triangle_count();
    out_mesh.resize_triangle(nfaces);
    std::size_t  f = 0;
    std::size_t max_vertex_id  = 0;
    for (tmesh_t::triangle_map_t::const_iterator it = mesh_.triangle_map().begin();
	 it != mesh_.triangle_map().end();  ++it) {
      const tmesh_t::triangle_t& t = it->first;
      out_mesh.set_triangle(f, std_triangle_mesh::tri_t(t[0],t[1],t[2]));
      for( std::size_t v =0; v <3; ++v ) {
	if ( t[v] > max_vertex_id  ) max_vertex_id = t[v];
      }
      f++;
    }
    std::cerr << "### Max vertex id " << max_vertex_id << std::endl;

    std::size_t nvert =  mesh_.vertex_count();
    std::cerr << "### Vertex count " << nvert << std::endl;

    out_mesh.resize_vertex(max_vertex_id + 1);
    for(std::size_t i = 0; i < max_vertex_id + 1; ++i) {
      const tmesh_t::vertex_data_t* v_ptr =  mesh_.vertex_data(i);
      if ( v_ptr != NULL ) { 
	out_mesh.set_position(i, v_ptr->position());
	tmesh_t::vertex_data_t::color_t c = v_ptr->color();
	out_mesh.set_color(i,std_triangle_mesh::color_t(c[0],c[1],c[2]));
	out_mesh.set_normal(i,v_ptr->normal());
	float q  = v_ptr->texel()[0];
	out_mesh.set_quality(i,q);
      } else {
	out_mesh.set_position(i, tmesh_t::vertex_data_t::point_t());
	out_mesh.set_color(i,std_triangle_mesh::color_t(0.0f,0.0f,0.0f));
	out_mesh.set_normal(i, tmesh_t::vertex_data_t::dual_vector_t(0.0f,0.0f,0.0f));	
      }
    }
    
    std::cerr<< "### Exported mesh with " << f << " triangles and " << nvert <<" valide vertices " <<  std::endl;
    
  }

  
}
