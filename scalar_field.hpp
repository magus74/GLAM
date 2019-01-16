#ifndef SCALAR_FIELD_HPP
#define SCALAR_FIELD_HPP

#include <vic/ply_io/std_triangle_mesh.hpp>
#include <sl/triangle_mesh.hpp>

namespace vic {

  // Vertex data decorated with scalar field phi and discrete differential geometry data ( gradient and laplacian)
  class vdata_t {
  public:
    typedef float        value_t;
    typedef sl::point3f  point_t;
    typedef sl::point2f  texel_t;
    typedef sl::vector3f vector_t;
    typedef vector_t     color_t;
    typedef sl::row_vector3f  dual_vector_t;
  protected:
    point_t   position_;
    texel_t   texel_;
    dual_vector_t normal_;
    vector_t  color_;
    std::size_t object_id_;
    
    value_t   phi_;
    vector_t  nabla_phi_;
    value_t   laplacian_phi_;
    
  public:
    vdata_t(const point_t& p = point_t(), const texel_t& t = texel_t()):position_(p), texel_(t) {
      phi_=0.0f;
      nabla_phi_=vector_t(0.0f,0.0f,0.0f);
      laplacian_phi_=0.0f;
      object_id_ = 0;
    }
    
  public: // Triangle Mesh interface

    inline point_t&  position()       { return position_; }
    inline const point_t&  position() const { return position_; }
    inline texel_t&  texel()       { return texel_; }
    inline const texel_t&  texel() const { return texel_; }
    inline const vector_t& color() const { return color_; }
    inline vector_t& color()  { return color_; }
    inline const dual_vector_t& normal() const { return normal_; }
    inline dual_vector_t& normal()  { return normal_; }
    inline value_t phi() const { return phi_; }
    inline value_t& phi()  { return phi_; }
    inline const vector_t& nabla_phi() const { return nabla_phi_; }
    inline vector_t& nabla_phi()  { return nabla_phi_; }
    inline value_t laplacian_phi() const { return laplacian_phi_; }
    inline value_t& laplacian_phi()  { return laplacian_phi_; }
    inline std::size_t object_id() const { return object_id_; }
    inline std::size_t& object_id()  { return object_id_; }


    inline vdata_t lerp(vdata_t& other, float t) const {
      vdata_t result = vdata_t();
      result.position() = position().lerp(other.position(), t);
      result.texel() = texel().lerp(other.texel(),t);
      result.color() = color().lerp(other.color(),t);
      result.normal() = normal().lerp(other.normal(),t);
      result.phi() = phi() + t*(other.phi()-phi());
      result.nabla_phi() = nabla_phi().lerp(other.nabla_phi(),t);
      result.laplacian_phi() = laplacian_phi() + t*(other.laplacian_phi()-laplacian_phi());
      result.object_id() = object_id();
      return result;
    }
    
  };
  
  struct edata_t {};

  struct tdata_t {};

  class scalar_field  {

  public:
    typedef sl::triangle_mesh<vdata_t, edata_t, tdata_t> tmesh_t; 
  
  public:
    void import_from(const vic::std_triangle_mesh& in_mesh);
    void export_to(vic::std_triangle_mesh& out_mesh);

    const tmesh_t& mesh()  const { return mesh_; }
    tmesh_t& mesh() { return mesh_; }

  protected:
    tmesh_t mesh_;

  };


}
#endif
