#ifndef INFLUENCE_MAP_BUILDER_HPP
#define INFLUENCE_MAP_BUILDER_HPP

#include <scalar_field.hpp>
#include <parser.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/oriented_box.hpp>

namespace vic {

  // FIXME --> more accurate transport theory needed
  class object_info {
  public:
    object_info():id_(0), alpha_max_(0.0),
      alpha_tot_(0.0), sigma_tot_(0.0), phi_max_(0.0), phi_tot_(0.0){}
  public:
    std::size_t id_;
    sl::aabox3f box_;
    double alpha_max_;
    double alpha_tot_;
    double sigma_tot_;
    double phi_max_;
    double phi_tot_;
  };


  class influence_map_builder {
  public:
    typedef vdata_t::value_t value_t;
    typedef vdata_t::point_t point_t;
    typedef vdata_t::vector_t vector_t;
    typedef vdata_t::color_t color_t;
    typedef scalar_field::tmesh_t mesh_t;

    typedef std::pair<sl::point3f, float> point_value_t;
    typedef std::vector< point_value_t> peak_vector_t;
    typedef std::pair<float,color_t> color_node_t;

    typedef enum { GAUSSIAN=0, LAMBERT} light_model_t;

  public:

    influence_map_builder(const parser& p): parser_(p) {
      update_max_weight_ = false;
      max_weight_ = 1000.0f;
      min_weight_ = 0.0f;
      max_laplacian_ = 1000.0f;
      max_gradient_ = 1000.0f;
      max_influence_distance_ = 0.1f;
      light_model_ = GAUSSIAN;
    }

    void set_parser(const parser& p) {
      parser_ = p;
    }

    void load_model(const std::string& model);
    void save_model(const std::string& model);

    void set_max_weight(float phi) {
      max_weight_  = phi;
    }
    
    void set_max_influence_distance(float phi) {
      max_influence_distance_ = phi;
    }

    void set_update_max_weight(bool b) {
      update_max_weight_ = b;
    }

    void set_light_model( light_model_t lm ) {
      light_model_ = lm;
    }

    void add_gradient_field(const std::vector<influence_map_builder::color_node_t>& cmap,
			    float sampling = 0.04f, // Vertex sampling
			    float threshold = 0.3f, // Length threshold
			    float max_length = 100.0f, // Max length
			    float max_radius = 10.0f); // Max radius

    void compute_color_map_field( const std::vector<color_node_t>& cmap,
				   bool laplacian_field = false, // Default scalar field
				   bool log_scale = false); // Default linear scale


    void colorize_objects(const std::vector<color_node_t>& cmap, float amax = 0.0f, float threshold = 0.1f);

    double compute_surface_flux();

    value_t max_weight() const { return max_weight_; }

    mesh_t& mesh() { return scalar_field_.mesh(); }

    void export_peak_positions(float cluster_threshold, 
			       std::vector< peak_vector_t> & peaks);

    const sl::aabox3f& bbox() const { return bbox_; }

    const std::vector<object_info>& obj_info() const { return object_info_; }
    
  protected:
    void create_mesh();
    void compute_scalar_field();
    void compute_gradient_field();
    void compute_laplacian();
    void colorize(const sl::interpolation_track<value_t, influence_map_builder::color_t>& color_map,
		  scalar_field::tmesh_t& mesh,
		  bool laplacian_field = false, // Default scalar field
		  bool log_scale = false);

    void label_objects();
    void compute_boxes();
    
    void export_obj(const std::string& model);

  protected:

    scalar_field scalar_field_;
    sl::aabox3f bbox_;
    parser parser_;

    value_t max_influence_distance_;
    value_t max_weight_;
    value_t max_laplacian_;
    value_t max_gradient_;
    value_t min_weight_;
    bool   update_max_weight_;
    light_model_t light_model_;

    std::vector<object_info> object_info_;
  };
}
#endif
