#ifndef VIC_STD_TRIANGLE_MESH_HPP
#define VIC_STD_TRIANGLE_MESH_HPP

#include <iostream>
#include <fstream>

#include <sl/fixed_size_point.hpp>
#include <vector>

namespace vic {

  class std_triangle_mesh {
  public:
    typedef float                         value_t;
    typedef sl::uint64_t 	           size_type_t;
    typedef sl::uint32_t                   uint_t;
    typedef sl::point3f                    point_t;
    typedef sl::point3f                    color_t;
    typedef sl::row_vector3f               dual_vector_t;
    typedef sl::fixed_size_point<3,uint_t> tri_t;

  public:
    std_triangle_mesh(uint_t nverts = 0,uint_t ntris = 0);
    ~std_triangle_mesh();

    void set_vertex(uint_t idx,const point_t& p=point_t(0.0,0.0,0.0),const color_t& c=color_t(1.0,1.0,1.0),const dual_vector_t& n=dual_vector_t(), float q=0.0f);
    void set_position(uint_t idx,const point_t& p);
    void set_color(uint_t idx,const color_t& c);
    void set_normal(uint_t idx,const dual_vector_t& n);
    void set_quality(uint_t idx, value_t q);
    void set_triangle(uint_t idx,const tri_t& t);


    point_t get_position(uint_t idx) const;
    color_t get_color(uint_t idx) const;
    dual_vector_t get_normal(uint_t idx) const;
    tri_t get_triangle(uint_t idx) const;
    value_t get_quality(uint_t idx) const;

    inline uint_t vertex_count() const;
    inline uint_t triangle_count() const;

    inline void resize_vertex(uint_t nverts);
    inline void resize_triangle(uint_t ntris);

    inline void push_vertex(const point_t& p=point_t(0.0,0.0,0.0),
			    const color_t& c=color_t(1.0,1.0,1.0),
			    const dual_vector_t& n=dual_vector_t(),
			    float q = 0.0f);

    inline void push_triangle(const tri_t& t);

  protected:
    bool valid_vertex_index(uint_t idx) const;
    bool valid_triangle_index(uint_t idx) const;
    
  protected:
    std::vector<point_t>        mesh_vertex_position_;
    std::vector<color_t>        mesh_vertex_color_;
    std::vector<dual_vector_t>  mesh_vertex_normal_;
    std::vector<value_t>        mesh_vertex_quality_;
    std::vector<tri_t>          mesh_triangles_;

  };

  inline std_triangle_mesh::uint_t std_triangle_mesh::vertex_count() const {
    return mesh_vertex_position_.size();
  }

  inline std_triangle_mesh::uint_t std_triangle_mesh::triangle_count() const {
    return mesh_triangles_.size();
  }

  inline void std_triangle_mesh::resize_triangle(uint_t ntris) {
    mesh_triangles_.resize(ntris);
  }

  inline void std_triangle_mesh::resize_vertex(uint_t nverts) {
    mesh_vertex_position_.resize(nverts);
    mesh_vertex_color_.resize(nverts);
    mesh_vertex_normal_.resize(nverts);
    mesh_vertex_quality_.resize(nverts);
  }

  inline void std_triangle_mesh::push_vertex(const point_t& p,
					     const color_t& c,
					     const dual_vector_t& n,
					     float q) {
    mesh_vertex_position_.push_back(p);
    mesh_vertex_color_.push_back(c);
    mesh_vertex_normal_.push_back(n);  
    mesh_vertex_quality_.push_back(q);
#if 0 
    std::cout << "push vertex - (" 
	      << p[0] << "," << p[1] << "," << p[2] << ") - ("
	      << c[0] << "," << c[1] << "," << c[2] << ") - ("
	      << n[0] << "," << n[1] << "," << n[2] << ")" << std::endl;
#endif
  }
  
  inline void std_triangle_mesh::push_triangle(const tri_t& t) {
    mesh_triangles_.push_back(t);
#if 0
    std::cout << "push triangle - (" 
	      << t[0] << "," << t[1] << "," << t[2] << ")" << std::endl;
#endif
  }


} //namespace vic

#endif

