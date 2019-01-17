#include <vic/ply_io/ooc_triangle_mesh.hpp>

namespace vic {

  ooc_triangle_mesh::ooc_triangle_mesh(uint_t nverts,uint_t ntris,std::string tmpfname,size_type_t cache_size) {
    mesh_vertex_position_ = new sl::external_array1<point_t>(tmpfname + "-positions.tmp","t",cache_size);
    mesh_vertex_color_    = new sl::external_array1<color_t>(tmpfname + "-colors.tmp","t",cache_size);
    mesh_vertex_normal_   = new sl::external_array1<dual_vector_t>(tmpfname + "-normals.tmp","t",cache_size);
    mesh_vertex_quality_   = new sl::external_array1<value_t>(tmpfname + "-qualities.tmp","t",cache_size);
    mesh_triangles_       = new sl::external_array1<tri_t>(tmpfname + "-triangles.tmp","t",cache_size);
    mesh_vertex_position_->resize(nverts);
    mesh_vertex_color_->resize(nverts);
    mesh_vertex_normal_->resize(nverts);
    mesh_vertex_quality_->resize(nverts);
    mesh_triangles_->resize(ntris);
  }

  void ooc_triangle_mesh::set_vertex(uint_t idx,const point_t& p,const color_t& c, const dual_vector_t& n, float q) {
    if(valid_vertex_index(idx)) {
      (*mesh_vertex_position_)[idx] = p;
      (*mesh_vertex_color_)[idx] = c;
      (*mesh_vertex_normal_)[idx] = n;
      (*mesh_vertex_quality_)[idx] = q;
    }
  }

  void ooc_triangle_mesh::set_position(uint_t idx,const point_t& p) {
     if(valid_vertex_index(idx)) {
      (*mesh_vertex_position_)[idx] = p;
     }
  }

  void ooc_triangle_mesh::set_color(uint_t idx,const color_t& c) {
    if(valid_vertex_index(idx)) {
      (*mesh_vertex_color_)[idx] = c;
    }
  }

  void ooc_triangle_mesh::set_normal(uint_t idx,const dual_vector_t& n) {
    if(valid_vertex_index(idx)) {
      (*mesh_vertex_normal_)[idx] = n;
    }
  }

  void ooc_triangle_mesh::set_quality(uint_t idx, value_t q) {
    if(valid_vertex_index(idx)) {
      (*mesh_vertex_quality_)[idx] = q;
    }
  }

  void ooc_triangle_mesh::set_triangle(uint_t idx,const tri_t& t) {
    if(valid_triangle_index(idx)) {
      (*mesh_triangles_)[idx] = t;
    }
  }


  ooc_triangle_mesh::point_t ooc_triangle_mesh::get_position(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return (*mesh_vertex_position_)[idx];
    }
    return point_t(0.0,0.0,0.0);
  }

  ooc_triangle_mesh::color_t ooc_triangle_mesh::get_color(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return (*mesh_vertex_color_)[idx];
    }
    return color_t(1.0,1.0,1.0);
  }

  ooc_triangle_mesh::dual_vector_t ooc_triangle_mesh::get_normal(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return (*mesh_vertex_normal_)[idx];
    }
    return dual_vector_t();
  }

  ooc_triangle_mesh::value_t ooc_triangle_mesh::get_quality(uint_t idx) const {
    if(valid_vertex_index(idx)) {
      return (*mesh_vertex_quality_)[idx];
    }
    return 0.0f;
  }

  ooc_triangle_mesh::tri_t ooc_triangle_mesh::get_triangle(uint_t idx) const {
    if(valid_triangle_index(idx)) {
      return (*mesh_triangles_)[idx];
    }
    return tri_t(0,0,0);
  }

  bool ooc_triangle_mesh::valid_vertex_index(uint_t idx) const {
    if(idx > (mesh_vertex_position_->size()-1)) {
      std::cerr << "ERROR - wrong vertex index" << std::endl;
      return false;
    }
    return true;
  }

  bool ooc_triangle_mesh::valid_triangle_index(uint_t idx) const {
    if(idx > (mesh_triangles_->size()-1)) {
      std::cerr << "ERROR - wrong triangle index" << std::endl;
      return false;
    }
    return true;
  }

  ooc_triangle_mesh::~ooc_triangle_mesh() {
    delete mesh_vertex_position_;
    mesh_vertex_position_ = 0;
    delete mesh_vertex_color_;
    mesh_vertex_color_ = 0;
    delete mesh_vertex_normal_;
    mesh_vertex_normal_ = 0;
    delete mesh_vertex_quality_;
    mesh_vertex_quality_ = 0;
    delete mesh_triangles_;
    mesh_triangles_ = 0;
    //std::cout << "tmp files deleted!!!" << std::endl;
  }
} //namespace vic
