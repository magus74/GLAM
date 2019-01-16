#ifndef VIC_CLUSTER_ANALYSIS_HPP
#define VIC_CLUSTER_ANALYSIS_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/fixed_size_vector.hpp>

namespace vic {

  class cluster_info {
  public:
    typedef float value_t;
    typedef sl::point3f point_t;
    typedef sl::vector3f vector_t;
    cluster_info() {
      pmin_ = point_t(1.0e9,1.0e9,1.0e9);
      pmax_ = point_t(-1.0e9,-1.0e9,-1.0e9);
      pmean_ = point_t();
      count_ = 0.0;
      volume_ = 0.0;
      vmin_ = 1.0e9;
      vmax_ = -1.0e9;
      vmean_ = 0.0;
    }
    
  public:
    point_t   pmin_;
    point_t   pmax_;
    point_t   pmean_;
    value_t   count_;
    value_t   volume_;
    value_t   vmin_;
    value_t   vmax_;
    value_t   vmean_;
  };
  
  
  class cluster_analysis {
  public:
    typedef float value_t;
    typedef sl::point3f point_t;
    typedef sl::vector3f vector_t;
    typedef std::pair<point_t,value_t> data_t;
    typedef std::vector<data_t> vector_data_t;
    typedef std::vector<point_t> point_data_t;
    
  public:
    cluster_analysis(const vector_data_t& d) {
      data_ = d;
    }
    
    void dbscan(float knn_maxradius = 0.03f, 
		std::size_t min_points = 100);
 
    const std::vector<cluster_info>& cluster_info_data() const { return cluster_info_; }

    void kmeans(std::size_t K, point_data_t& seeds);

  protected:
    float kmeans_assign(const point_data_t& centers);
    void kmeans_centroid(point_data_t& centers);
    void fill_cluster_barycenters(std::size_t k, bool knn = true);
    void compute_cluster_info(std::vector<cluster_info>& ci, std::size_t k, bool knn);

  protected:    
    std::vector<int> id_;
    vector_data_t data_;
    point_data_t cluster_centers_;
    std::vector<cluster_info> cluster_info_;
  };
  
}

#endif
