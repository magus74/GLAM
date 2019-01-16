#include <parser.hpp>
#include <vic/ply_io/ooc_ply_io.hpp>

namespace vic {
  
  void parser::create_kdtree(const std::string& srcfile) {
    ooc_triangle_mesh * mesh = ooc_ply_io::from_file(srcfile.c_str());
    // Compute bounding box
    for (std::size_t i = 0; i < mesh->vertex_count(); ++i) {
      vic::ooc_triangle_mesh::point_t p = mesh->get_position(i);
      kdtree_.insert(point_t(p[0],p[1],p[2]));
    }
    std::cerr << " Created kdtree with " << mesh->vertex_count() << " vertices" << std::endl;
  }
    
}
