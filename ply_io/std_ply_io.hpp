#ifndef VIC_STD_PLY_IO_HPP
#define VIC_STD_PLY_IO_HPP

#include <vic/ply_io/std_triangle_mesh.hpp>
#include <sl/rigid_body_map.hpp>

namespace vic {

  class std_ply_io {
  public:
    typedef sl::rigid_body_map<3,float> rigid_body_map_t;
    typedef std_triangle_mesh::point_t  point_t;

  public:     
    static void from_file(const char* str, std_triangle_mesh& tm, float scale = 1.0f);
    static void to_file(const std_triangle_mesh& tm, const char* str);

    static void apply_transform(const char* str_dest,
				const char* str_src,
				rigid_body_map_t& R,
				rigid_body_map_t& T);
   
  };

} // namespace vic

#endif
