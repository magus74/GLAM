#include <influence_map_builder.hpp>
#include <vic/ply_io/std_ply_io.hpp>
#include <sl/clock.hpp>
#include <sl/random.hpp>
#include <list>
#include <fstream>

namespace vic {

	static const double M_PI = 3.1415;

  static inline float gaussian(float x2, float xmax2) {
    float result=0.0f;
    if (x2<xmax2) {
      float xx = x2/xmax2;
      result = 1.0f-xx;
      result *= result;
    }
    return result;
  }

  static double surface_area(const scalar_field::tmesh_t&mesh,
			     const scalar_field::tmesh_t::triangle_t& tri) {
    const vdata_t* v_i = mesh.vertex_data(tri[0]);
    const vdata_t* v_j = mesh.vertex_data(tri[1]);
    const vdata_t* v_k = mesh.vertex_data(tri[2]);
    
    vdata_t::vector_t e_ij = v_j->position() - v_i->position();
    vdata_t::vector_t e_ik = v_k->position() - v_i->position();
    return 0.5*(e_ij.cross(e_ik).two_norm());
  }

  static double surface_area(const scalar_field::tmesh_t&mesh,
			     const scalar_field::tmesh_t::small_triangle_set_t& tri_set) {

    double result = 0.0;
    for (scalar_field::tmesh_t::small_triangle_set_t::const_iterator t_it = tri_set.begin();
             t_it != tri_set.end();
             ++t_it) {
      result += surface_area( mesh, *t_it);
    }
    return result;
  }

  static double surface_area(const scalar_field::tmesh_t&mesh ) {
    typedef scalar_field::tmesh_t::const_triangle_iterator_t const_tri_it;
    double result = 0.0;
    for (const_tri_it  t_it = mesh.triangle_begin(); t_it != mesh.triangle_end(); ++ t_it) {
      result += surface_area(mesh, t_it->first);
    } 
    return result;    
  }

  static void push_disk(scalar_field::tmesh_t&mesh, 
			const vdata_t::point_t&v, 
			const vdata_t::vector_t&d, 
			const vdata_t::color_t&c, 
			float r, 
			std::size_t tess) {

    typedef  scalar_field::tmesh_t mesh_t;
    
    vdata_t::vector_t x_hat = vdata_t::vector_t(0.0f,1.0f,0.0f).cross(-d).ok_normalized();
    vdata_t::vector_t y_hat = x_hat.cross(d).ok_normalized();

    mesh_t::vertex_t start_id = mesh.vertex_count(); 
    for(std::size_t i = 1; i<= tess; ++i) {
      mesh.insert_triangle ( mesh_t::triangle_t(start_id, start_id + i, start_id + (1+i)%tess));
      float theta = 2.0f * float(i-1)/ float(tess)  * M_PI; 
      mesh_t::vertex_data_t* v_ptr =  mesh.vertex_data(start_id + i);
      v_ptr->position() = v + r*cosf(theta)*x_hat + r*sinf(theta)*y_hat;
      v_ptr->color() = c;
      v_ptr->normal() =  vdata_t::dual_vector_t(d[0],d[1],d[2]);
    }

    mesh_t::vertex_data_t* v_ptr =  mesh.vertex_data(start_id);
    v_ptr->position() = v;
    v_ptr->color() = c;
    v_ptr->normal() =   vdata_t::dual_vector_t(d[0],d[1],d[2]);

  }

  static void push_arrow_glyph(scalar_field::tmesh_t& mesh,
			       const vdata_t::point_t& p,
			       const vdata_t::vector_t& d,
			       const vdata_t::color_t& c,
			       float /*l*/, float r) {
    
    const std::size_t tess = 10;
    push_disk(mesh, p, d, c, 0.5f*r, tess);
#if 0
    push_cylinder(mesh, p, d, 0.5f*r, 0.8f*l, tess);
    vdata_t::point_t v_t = p + 0.8f * l * d;
    push_ring(mesh, v_t, d, 0.5f*r, r, tess);
    push_cone(mesh, v_t, d, 0.2f*l, r, tess);
#endif
  }


  typedef sl::axis_aligned_box_builder<3, float> abuilder3f_t;

  static void compute_axis_aligned_box( const scalar_field& sf, sl::aabox3f& b, std::size_t id = 0) {
    const scalar_field::tmesh_t& mesh = sf.mesh();
    abuilder3f_t builder;
    builder.begin_model();
    for( scalar_field::tmesh_t::const_vertex_iterator_t  it = mesh.vertex_begin(); it != mesh.vertex_end(); ++it) {
      const scalar_field::tmesh_t::vertex_data_t* v_ptr =  it->second.data();
      if ( v_ptr->object_id() == id ) builder.put_point(v_ptr->position());
    }
    builder.end_model();
    b = builder.last_bounding_volume();
  }


  void influence_map_builder::load_model(const std::string& model) {
    
    std_triangle_mesh  mesh;
    std_ply_io::from_file(model.c_str(), mesh);

    scalar_field_.import_from(mesh);
    compute_axis_aligned_box(scalar_field_, bbox_, 0);

    std::cerr<< "# Bbox min: "<< bbox_[0][0] <<" ,"<< bbox_[0][1] <<" ,"<< bbox_[0][2] << std::endl;
    std::cerr<< "# Bbox max: "<< bbox_[1][0] <<" ,"<< bbox_[1][1] <<" ,"<< bbox_[1][2] << std::endl;

    label_objects();
    compute_boxes();
    
  }

  void influence_map_builder::export_obj(const std::string& model) {
    std::ofstream out_file(model.c_str());
    // Header with material reference
    if (model.find("_glam") != std::string::npos ) {
      out_file <<"mtllib glam.mtl"<<std::endl;
    } else if (model.find("_peak") != std::string::npos ) {
      out_file <<"mtllib peak.mtl"<<std::endl;
    } else if (model.find("_obj") != std::string::npos ) {
      out_file <<"mtllib obj.mtl"<<std::endl;
    } else {
      std::cerr << "### Warning material file not recognized "<< std::endl;
      out_file <<"mtllib default.mtl"<<std::endl;
    }
    out_file <<"usemtl mat"<<std::endl;
    typedef scalar_field::tmesh_t mesh_t;
    mesh_t& mesh = scalar_field_.mesh();
    std::map<std::size_t,std::size_t> index_map;
    std::size_t idx = 1;
    for( mesh_t::const_vertex_iterator_t v_it = mesh.vertex_begin(); v_it != mesh.vertex_end(); ++v_it) {
      std::size_t v_id = v_it->first;
      const vdata_t* v_i = mesh.vertex_data(v_id);
      vdata_t::point_t p = v_i->position();
      vdata_t::texel_t t = v_i->texel();
      vdata_t::dual_vector_t n = v_i->normal(); 
      out_file << "v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
      out_file << "vt "<<t[0]<<" "<<t[1]<<std::endl;
      out_file << "vn "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
      index_map.insert( std::make_pair(v_id,idx));
      idx++;
    }
    for( mesh_t::const_triangle_iterator_t f_it = mesh.triangle_begin(); f_it != mesh.triangle_end(); ++f_it) {
      const mesh_t::triangle_t& tr = f_it->first;
      if ( index_map.find(tr[0]) == index_map.end() ||
	   index_map.find(tr[1]) == index_map.end() ||
	   index_map.find(tr[2]) == index_map.end() ) {
	continue;
      }
      std::size_t i = index_map[ tr[0] ];
      std::size_t j = index_map[ tr[1] ];
      std::size_t k = index_map[ tr[2] ];
      out_file <<"f "<<i<<"/"<<i<<"/"<<i
	       <<" "<<j<<"/"<<j<<"/"<<j
	       <<" "<<k<<"/"<<k<<"/"<<k<< std::endl;
    }
    out_file.close();
  }

  void influence_map_builder::save_model(const std::string& model) {

    bool ply_output = model.find(".ply") != std::string::npos;
    if ( ply_output ) {
      std_triangle_mesh mesh;
      scalar_field_.export_to(mesh);
      vic::std_ply_io::to_file(mesh, model.c_str());
    } else {
      // Assuming wavefront obj output
      if ( model.find(".obj") == std::string::npos ) {
	std::cerr <<" ### Error: wrong output file name  in " << model << std::endl;
	return;
      } 
      export_obj(model);
    }
    std::cerr <<"# Written file " << model << std::endl;

  }

  void influence_map_builder::compute_boxes() {
    for (std::size_t i = 0; i < object_info_.size(); ++i) {
      compute_axis_aligned_box( scalar_field_, object_info_[i].box_, i );
    }
    std::cerr <<"# Computed boxes for " << object_info_.size() << " objects " << std::endl;
  }

  void influence_map_builder::label_objects() {

    typedef scalar_field::tmesh_t mesh_t;
    mesh_t& mesh = scalar_field_.mesh();
    object_info_.clear();
   
    std::set<mesh_t::vertex_t> labelled;
    std::size_t current_object = 0;
    for( mesh_t::vertex_iterator_t v_it = mesh.vertex_begin(); v_it != mesh.vertex_end(); ++v_it) {
      mesh_t::vertex_t vid_i = v_it->first;
      if ( labelled.find(vid_i) == labelled.end() ) {
	std::list<mesh_t::vertex_t> to_check;
	object_info info;
	info.id_ = current_object;
	to_check.push_back(vid_i);
	//	std::cerr<<" Obj " << current_object << " first push " << vid_i << std::endl;
	while ( to_check.size() ) {
	  mesh_t::vertex_t cur_ver = to_check.front();
	  // label current vertex
	  vdata_t* v_it = mesh.vertex_data(cur_ver);
	  v_it->object_id() = current_object;
	  labelled.insert(cur_ver);
	  to_check.pop_front();
	  //std::cerr<<" Obj " << current_object << " pop " << cur_ver << std::endl;

	  // add connected vertices
	  const mesh_t::small_edge_set_t* v_edges = mesh.vertex_edges(cur_ver);
	  if ( mesh.vertex_edge_count(cur_ver) < 2 ) {
	    // std::cerr << " Error in vertex " << cur_ver << std::endl;
	  }
	  for (mesh_t::small_edge_set_t::const_iterator e_it = v_edges->begin(); e_it != v_edges->end(); ++e_it) {
	    mesh_t::edge_t e = *e_it;
	    mesh_t::vertex_t vid_j = ( e[0] != cur_ver ? e[0] : e[1]);
	    //std::cerr<<" Obj " << current_object << " check " << vid_j << std::endl;
	    if (labelled.find(vid_j) == labelled.end()) {
	      to_check.push_back(vid_j);
	      vdata_t* vj_it = mesh.vertex_data(vid_j);
	      vj_it->object_id() = current_object;
	      labelled.insert(vid_j);
	      //std::cerr<<" Obj " << current_object << " push " << vid_j << std::endl;
	    }
	  }
	}
	object_info_.push_back(info);
	current_object++;
      }
    }
    std::cerr<<" Labelled model with " << current_object << " objects " << std::endl; 
  }



  void influence_map_builder::export_peak_positions(float cluster_threshold, 
						    std::vector< influence_map_builder::peak_vector_t  >& peaks) {
    typedef scalar_field::tmesh_t mesh_t;
    const mesh_t& mesh = scalar_field_.mesh();
    peaks.resize( object_info_.size());

    for( mesh_t::const_vertex_iterator_t v_it = mesh.vertex_begin(); v_it != mesh.vertex_end(); ++v_it) {
      mesh_t::vertex_t vid_i = v_it->first;
      const vdata_t* v_i = mesh.vertex_data(vid_i);
      value_t phi = v_i->phi();
      if ( phi >= cluster_threshold * max_weight_ ) {
	if ( v_i->object_id() < peaks.size() ) {
	  peaks[v_i->object_id()].push_back( std::make_pair(v_i->position(),phi) );
	}
      } 
    }
    std::cerr <<" Computed peaks for " << peaks.size() << " objects " << std::endl;
  }

  double influence_map_builder::compute_surface_flux() {
    typedef scalar_field::tmesh_t mesh_t;
    typedef mesh_t::const_triangle_iterator_t const_tri_it;
    typedef mesh_t::triangle_t triangle_t;

    const mesh_t& mesh = scalar_field_.mesh();
    double flux_tot = 0.0;
    double flux_max = 0.0;
#if 0
    double phi_tot  = 0.0;
    for( mesh_t::const_vertex_iterator_t v_it = mesh.vertex_begin(); v_it != mesh.vertex_end(); ++v_it) {
      mesh_t::vertex_t vid_i = v_it->first;
      const vdata_t* v_i = mesh.vertex_data(vid_i);
      value_t phi = v_i->phi();
      phi_tot += phi;
      const mesh_t::small_triangle_set_t& v_tri  =  v_it->second.triangle_set();
      double flux_local = phi * surface_area(mesh, v_tri);
      if ( flux_max  < flux_local ) flux_max = flux_local;
      flux_tot += flux_local;
    }
    double surface_tot = surface_area(mesh);
#else
    double surface_tot = 0.0;
    for (const_tri_it  t_it = mesh.triangle_begin(); t_it != mesh.triangle_end(); ++ t_it) {
      
      const triangle_t& tri = t_it->first;
      const vdata_t* v_i = mesh.vertex_data(tri[0]);
      const vdata_t* v_j = mesh.vertex_data(tri[1]);
      const vdata_t* v_k = mesh.vertex_data(tri[2]);

      double phi_tri = 1.0/3.0 *(v_i->phi()+v_j->phi()+v_k->phi());
      vdata_t::vector_t e_ij = v_j->position() - v_i->position();
      vdata_t::vector_t e_ik = v_k->position() - v_i->position();
      double dsigma = 0.5*(e_ij.cross(e_ik).two_norm());

      double flux_local = phi_tri * dsigma;
      if ( flux_max < flux_local ) flux_max = flux_local;
      surface_tot += dsigma;
      flux_tot += flux_local;

      std::size_t cur_id = v_i->object_id();
      if ( cur_id < object_info_.size() ) {
	object_info& info  = object_info_[v_i->object_id()];
	if ( info.id_ != cur_id ) {
	  std::cerr<< "### Error: mismatch between ids and vector order: " << info.id_ << " vs " << cur_id; 
	} else { 
	  if ( phi_tri > info.alpha_max_ ) info.alpha_max_ = phi_tri;
	  if ( flux_local > info.phi_max_ ) info.phi_max_ = flux_local;
	  info.alpha_tot_ += phi_tri;
	  info.sigma_tot_ += dsigma;
	  info.phi_tot_ += flux_local;
	}
      } else {
	  std::cerr<< "### Error: object id not correctly labelled: " << cur_id; 
      }

    }
#endif

    std::cerr <<" ##  Surface total == " << surface_tot << std::endl;
    std::cerr <<" ##  Flux maximum  == " << flux_max << std::endl;
    std::cerr <<" ##  Flux total  == " << flux_tot << std::endl;

    return flux_tot / surface_tot; 
  }
  
  void influence_map_builder::compute_gradient_field() {
    typedef scalar_field::tmesh_t mesh_t;
    mesh_t& mesh = scalar_field_.mesh();
    float max_nabla = 0.0f;
    for( mesh_t::vertex_iterator_t v_it = mesh.vertex_begin(); v_it != mesh.vertex_end(); ++v_it) {
      mesh_t::vertex_t vid_i = v_it->first;
      vdata_t* v_i = mesh.vertex_data(vid_i);
      mesh_t::small_edge_set_t v_edges = v_it->second.edge_set();
      for (mesh_t::small_edge_set_t::const_iterator e_it = v_edges.begin(); e_it != v_edges.end(); ++e_it) {
	mesh_t::edge_t e = *e_it;
	mesh_t::vertex_t vid_j = ( e[0] != vid_i ? e[0] : e[1]);
	vdata_t* v_j = mesh.vertex_data(vid_j);
	vdata_t::vector_t dp = v_i->position() - v_j->position();
	if ( dp.two_norm() > 1.0e-6) {
	  v_i->nabla_phi() += vdata_t::vector_t(1.0f/(dp[0]+1.0e-6),1.0f/(dp[1]+1.0e-6),1.0f/(dp[2]+1.0e-6));
	}
      }
      v_i->nabla_phi() *= v_i->phi();
      float nabla = v_i->nabla_phi().two_norm();
      if (nabla > max_nabla) {
	max_nabla = nabla;
      }
    }
    max_gradient_ = max_nabla;
    std::cerr<<"#### Computed gradient with max value " << max_gradient_ <<  std::endl;
  }

  void influence_map_builder::compute_laplacian() {
    typedef scalar_field::tmesh_t mesh_t;
    mesh_t& mesh = scalar_field_.mesh();
    for( mesh_t::edge_iterator_t e_it = mesh.edge_begin(); e_it != mesh.edge_end(); ++e_it) {
      mesh_t::edge_t e = e_it->first;
      vdata_t* v_i = mesh.vertex_data(e[0]);
      vdata_t* v_j = mesh.vertex_data(e[1]);
      mesh_t::edge_attribute_t e_attr = e_it->second;
      mesh_t::small_triangle_set_t tri_edges = e_attr.triangle_set();
      for (mesh_t::small_triangle_set_t::const_iterator t_it = tri_edges.begin(); t_it != tri_edges.end(); ++t_it) {
	mesh_t::triangle_t t = *t_it;
	bool found = false;
	vdata_t* v_k;
	for(std::size_t k = 0; k< 3; ++k) {
	  found = ( (t[k] != e[0]) && (t[k] != e[1]) );
	  if (found) {
	    v_k = mesh.vertex_data(t[k]);
	    break;
	  }
	}
	if (!found) {
	  std::cerr<< "(EEE): Problem with mesh " << std::endl;
	  std::cerr<< " Triangle " << t[0] <<", "<< t[1]<<", "<<t[2]<<std::endl;
	  exit(1);
	} 
	float alpha_ij = fabs( (v_i->position() - v_k->position()).angle(v_j->position()-v_k->position()));
	if ( alpha_ij > 1.0e-10 ) {
	  float w_ij = 0.5f * cosf(alpha_ij)/sinf(alpha_ij) * (v_i->phi() - v_j->phi());
	  v_i->laplacian_phi() += w_ij;
	  v_j->laplacian_phi() -= w_ij;
	}
      }
    }
    // Compute max laplacian
    std::size_t nvert =  mesh.vertex_count();
    float max_lap = 0.0f;
    for( std::size_t  i = 0; i < nvert; ++i) {
      scalar_field::tmesh_t::vertex_data_t* v_ptr =  mesh.vertex_data(i);
      if ( v_ptr->laplacian_phi() > max_lap ) {
	max_lap = v_ptr->laplacian_phi();
      }
    }
    max_laplacian_ = max_lap;
    std::cerr<<"#### Computed laplacian with max value == " << max_laplacian_<< std::endl;
  }
 
  void influence_map_builder::compute_scalar_field() {

    std::cerr<<"### Start computing scalar field" << std::endl;
    sl::real_time_clock ck; 
    scalar_field::tmesh_t& mesh = scalar_field_.mesh();
    std::size_t nvert =  mesh.vertex_count();
    float max_phi = 0.0f;
    std::size_t i = 0;
    double phi_total = 0.0f;
    for( scalar_field::tmesh_t::vertex_iterator_t  it = mesh.vertex_begin(); it != mesh.vertex_end(); ++it) {
      scalar_field::tmesh_t::vertex_data_t* v_ptr =  it->second.data();
      sl::point3f p = v_ptr->position();
      vdata_t::dual_vector_t n = v_ptr->normal();
      std::vector<parser::point_t> knn;
      parser_.nearest_neighbors( knn, p, max_influence_distance_);
      // Integration with light decay model
      for (std::size_t j = 0; j < knn.size(); ++j ) {
	sl::vector3f d_j = knn[j] - p;
	float tau = d_j.ok_normalized().dot(sl::vector3f(n[0],n[1],n[2]));
	float d2 =  d_j.two_norm_squared();
	if ( tau >= 0.0f ) {
	  if ( light_model_ == GAUSSIAN ) {
#if 1
	    v_ptr->phi() += gaussian( d2, 
				    max_influence_distance_*max_influence_distance_);
#else
	    v_ptr->phi() += gaussian( sqrtf(d2), 
				    max_influence_distance_);
#endif
	  } else { // LAMBERT
	    v_ptr->phi() += tau * (1.0f - d2 );
	  }
	}
      }
      if ( update_max_weight_) {
	if ( v_ptr->phi() > max_phi ) {
	  max_phi = v_ptr->phi();
	}
      }
      if (++i > 0 &&  i%10000  == 0) {
	float eta_sec = float(ck.elapsed().as_seconds())*float(nvert/float(i) -1.0f);
	int eta_min = int(eta_sec/60.0f);
	int eta_hour = int(eta_min/60.0f);
	std::cout<<"# V: "<< i << "/" << nvert <<", eta: " << eta_hour 
		 <<"h "<<eta_min%60<<"m "<<int(eta_sec)%60<<"s"<<  std::endl;
      }
      phi_total += double(v_ptr->phi());
    }
    if ( update_max_weight_) {
      max_weight_ = max_phi;
      std::cerr<<"### Updated max weight == " << max_weight_ <<  std::endl;
    }
    std::cerr<<"### Finished computing scalar field" << std::endl;
    std::cerr<<"### Total scalar value   == " << phi_total << std::endl;
  }

  // Add gradient field per vertex ==> for the moment we represent it as glyph
  void influence_map_builder::add_gradient_field(const std::vector<influence_map_builder::color_node_t>& cmap,
					     float sampling, // Vertex sampling
					     float threshold, // Length threshold
					     float max_length, // Max length
					     float max_radius) {
    compute_gradient_field();
    scalar_field::tmesh_t& mesh = scalar_field_.mesh();
    std::size_t nvert =  mesh.vertex_count();
    
    // Create interpolated color map
    sl::interpolation_track<value_t,color_t> color_map;
    color_map.set_continuity(0);
    for(std::size_t i=0; i< cmap.size(); ++i) {
      color_map.insert( cmap[i]);
    }     
    
    std::size_t ratio = std::size_t(1.0f/threshold);
    for( std::size_t  i = 0; i < nvert; i+= ratio) {
      scalar_field::tmesh_t::vertex_data_t* v_ptr =  mesh.vertex_data(i);
      vdata_t::vector_t nabla = v_ptr->nabla_phi();
      vdata_t::value_t nabla_norm = nabla.two_norm();
      float csi = log(nabla_norm +1)/log(max_gradient_+1);
      if (csi > threshold) {
	color_t c = color_map.value_at(csi);
	vdata_t::value_t l = csi * max_length;
	vdata_t::value_t r = csi * max_radius;
	push_arrow_glyph(mesh, v_ptr->position(), (1.0f / nabla_norm)*nabla, c, l, r);
      }
    } 
  }
  
					     
  void influence_map_builder::compute_color_map_field(const std::vector<influence_map_builder::color_node_t>& cmap,
						  bool laplacian_field, // Default scalar field
						  bool log_scale) {

    // Use log parser for computing influence map according to the views
    compute_scalar_field();

    if (laplacian_field) {
      compute_laplacian();
    }

    // Create interpolated color map
    sl::interpolation_track<value_t,color_t> color_map;
    color_map.set_continuity(0);
    for(std::size_t i=0; i< cmap.size(); ++i) {
      color_map.insert( cmap[i]);
    }

    colorize(color_map, scalar_field_.mesh(), 
	     laplacian_field,  log_scale);
  }

  void influence_map_builder::colorize_objects(const std::vector<influence_map_builder::color_node_t>& cmap, float amax, float threshold) {

    // Computed expected absorptions
    std::vector<float> alpha(object_info_.size());

    float amin = threshold * max_weight_;
    std::set<std::size_t> objects_to_filter;
    for (std::size_t i =0; i < object_info_.size(); ++i ) {
      if ( object_info_[i].alpha_max_ < amin ) {
	objects_to_filter.insert(i);
      } 
      alpha[i] =  object_info_[i].phi_tot_ /  object_info_[i].sigma_tot_;
      amax = sl::max( amax, alpha[i] );
    }
    std::cerr<<"## Objects below threshold: " << objects_to_filter.size() << std::endl;
    std::cerr <<" # Using maximum expected absorption " << amax << std::endl; 
   
    // Create interpolated color map
    sl::interpolation_track<value_t,color_t> color_map;
    color_map.set_continuity(0);
    for(std::size_t i=0; i< cmap.size(); ++i) {
      color_map.insert( cmap[i]);
    }     

    scalar_field::tmesh_t& mesh = scalar_field_.mesh();
    std::set<scalar_field::tmesh_t::vertex_t> verts_to_remove;

    //std::size_t nvert =  mesh.vertex_count();   
    for( scalar_field::tmesh_t::vertex_iterator_t  it = mesh.vertex_begin(); it != mesh.vertex_end(); ++it) {
      scalar_field::tmesh_t::vertex_data_t* v_ptr =  it->second.data();
      // v_ptr->color() =  color_map.value_at(t);
      std::size_t id = v_ptr->object_id();
      if ( objects_to_filter.find(id) != objects_to_filter.end() ) {
	verts_to_remove.insert( it->first );
      } else {
	float t = (amax > 0.0f ? alpha[id]/amax: 0.0f); 
	v_ptr->texel() = vdata_t::texel_t(sl::median(0.001f,t,1.0f),0.5f);
	v_ptr->color() = color_map.value_at(sl::median(0.0f,t,1.0f));
      }
    }

    for( std::set<scalar_field::tmesh_t::vertex_t>::iterator it = verts_to_remove.begin(); it != verts_to_remove.end(); ++it ) {
      if (mesh.has_vertex(*it)) { 
	mesh.erase_vertex(*it);
      }
    }

    std::cerr <<"### Colorized objects and filtered mesh" << std::endl;
  }

  void influence_map_builder::colorize(const sl::interpolation_track<value_t, influence_map_builder::color_t>& color_map,
				   scalar_field::tmesh_t& mesh,
				   bool laplacian_field,
				   bool log_scale) {

    //    std::size_t nvert =  mesh.vertex_count();   
    std::vector<color_t> obj_col;
    sl::random::unit_normal<float> rnd;
    for (std::size_t i = 0; i < object_info_.size(); ++i) {
      obj_col.push_back(color_t(rnd.value(), rnd.value(),rnd.value()));  
    }

    for( scalar_field::tmesh_t::vertex_iterator_t  it = mesh.vertex_begin(); it != mesh.vertex_end(); ++it) {
      scalar_field::tmesh_t::vertex_data_t* v_ptr =  it->second.data();
      float t =  v_ptr->phi()/max_weight_; 
      if (log_scale) {
	t = log(1.0f+ v_ptr->phi())/log(1.0f+max_weight_);
      }
      // v_ptr->color() =  color_map.value_at(t);

#if 1
      t = sl::median(0.001f,t,1.0f);
      v_ptr->color() = color_map.value_at(t);
      v_ptr->texel() = scalar_field::tmesh_t::vertex_data_t::texel_t(t,0.5f);
#else
      std::size_t oid = v_ptr->object_id();
      v_ptr->color() = 0.8f *color_map.value_at(t) + 0.2f *  obj_col[oid]; 
#endif
      if (laplacian_field) {
	t = v_ptr->laplacian_phi()/max_laplacian_;
	if (log_scale) {
	  t =  log(1.0f+ v_ptr->laplacian_phi())/log(1.0f+max_laplacian_);
	  t =  (t>0.7f?t:0.0f);
	}
	v_ptr->color() *= (1.0f - t);
      }
    }
    std::cerr <<"### Scalar function color mapped " << std::endl;
  }

}
