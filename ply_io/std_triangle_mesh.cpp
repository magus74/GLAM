#include <vic/ply_io/std_triangle_mesh.hpp>

namespace vic {

  std_triangle_mesh::std_triangle_mesh(uint_t nverts,uint_t ntris) {
    mesh_vertex_position_.resize(nverts);
    mesh_vertex_color_.resize(nverts);
    mesh_vertex_normal_.resize(nverts);
    mesh_vertex_quality_.resize(nverts);
    mesh_triangles_.resize(ntris);
  }

  void std_triangle_mesh::set_vertex(uint_t idx,const point_t& p,const color_t& c, const dual_vector_t& n, float q) {
    if(valid_vertex_index(idx)) {
      mesh_vertex_position_[idx] = p;
      mesh_vertex_color_[idx] = c;
      mesh_vertex_normal_[idx] = n;
      mesh_vertex_quality_[idx] = q;
    }
  }

  void std_triangle_mesh::set_position(uint_t idx,const point_t& p) {
     if(valid_vertex_index(idx)) {
      mesh_vertex_position_[idx] = p;
     }
  }

  void std_triangle_mesh::set_color(uint_t idx,const color_t& c) {
    if(valid_vertex_index(idx)) {
      mesh_vertex_color_[idx] = c;
    }
  }

  void std_triangle_mesh::set_normal(uint_t idx,const dual_vector_t& n) {
    if(valid_vertex_index(idx)) {
      mesh_vertex_normal_[idx] = n;
    }
  }

  void std_triangle_mesh::set_quality(uint_t idx, value_t q) {
    if(valid_vertex_index(idx)) {
      mesh_vertex_quality_[idx] = q;
    }
  }

  void std_triangle_mesh::set_triangle(uint_t idx,const tri_t& t) {
    if(valid_triangle_index(idx)) {
      mesh_triangles_[idx] = t;
    }
  }


  std_triangle_mesh::point_t std_triangle_mesh::get_position(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return mesh_vertex_position_[idx];
    }
    return point_t(0.0,0.0,0.0);
  }

  std_triangle_mesh::color_t std_triangle_mesh::get_color(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return mesh_vertex_color_[idx];
    }
    return color_t(1.0,1.0,1.0);
  }

  std_triangle_mesh::dual_vector_t std_triangle_mesh::get_normal(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return mesh_vertex_normal_[idx];
    }
    return dual_vector_t();
  }

  std_triangle_mesh::value_t std_triangle_mesh::get_quality(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return mesh_vertex_quality_[idx];
    }
    return 0.0f;
  }

  std_triangle_mesh::tri_t std_triangle_mesh::get_triangle(uint_t idx) const {
    if(valid_triangle_index(idx)) {
      return mesh_triangles_[idx];
    }
    return tri_t(0,0,0);
  }

  bool std_triangle_mesh::valid_vertex_index(uint_t idx) const {
    if(idx > (mesh_vertex_position_.size()-1)) {
      std::cerr << "ERROR - wrong vertex index" << std::endl;
      return false;
    }
    return true;
  }

  bool std_triangle_mesh::valid_triangle_index(uint_t idx) const {
    if(idx > (mesh_triangles_.size()-1)) {
      std::cerr << "ERROR - wrong triangle index" << std::endl;
      return false;
    }
    return true;
  }

  std_triangle_mesh::~std_triangle_mesh() {
    //std::cout << "tmp files deleted!!!" << std::endl;
  }
} //namespace vic
