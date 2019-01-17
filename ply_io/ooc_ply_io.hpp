#ifndef VIC_OOC_PLY_IO_HPP
#define VIC_OOC_PLY_IO_HPP

#include <vic/ply_io/ooc_triangle_mesh.hpp>
#include <sl/rigid_body_map.hpp>

namespace vic {

  class ooc_ply_io {
  public:
    typedef sl::rigid_body_map<3,float> rigid_body_map_t;
    typedef ooc_triangle_mesh::point_t  point_t;

  public:     
    static ooc_triangle_mesh* from_file(const char* str, float scale = 1.0f);
    static void to_file(ooc_triangle_mesh* tm, const char* str);

    static void apply_transform(const char* str_dest,
				const char* str_src,
				rigid_body_map_t& R,
				rigid_body_map_t& T);
   
  };

} // namespace vic

#endif
