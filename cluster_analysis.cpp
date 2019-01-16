#include <cluster_analysis.hpp>
#include <sl/fixed_size_vector.hpp>
#include <sl/random.hpp>
#include <sl/kdtree.hpp>
#include <fstream>

namespace vic {


  /*********************************************
Cluster view positions by using algorithm
DBSCAN: A density-based algorithm for discovering clusters in large spatial databases with noise. 
Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96)
Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996)
  **********************************************/
  void cluster_analysis::dbscan(float knn_maxradius, 
				std::size_t min_points) {

    typedef std::map<point_t, std::pair<std::size_t, int> > cluster_t;    
    std::size_t data_count = data_.size();
    id_.resize(data_count);
    sl::kdtree<3,float,point_t> kdt;
    cluster_t db_map;
    // Create kd-tree 
    for (std::size_t s = 0; s < data_count; ++s) {
      point_t csi = data_[s].first;
      db_map.insert ( std::make_pair( csi, std::make_pair(s, -1)) ); // Mark not visited
      kdt.insert(csi);
    }
  
    std::size_t cluster_id = 0;
    std::size_t visited_points = 0;
    for (cluster_t::iterator it = db_map.begin(); it != db_map.end(); ++it) {
      if ( (it->second).second < 0 ) { // not visited
	point_data_t near_points;
	kdt.k_nearest_neighbors_in(near_points, it->first, 1000, knn_maxradius); 
	if (near_points.size() < min_points ) {
	  (it->second).second = 0; //  View point marked as noise
	  ++visited_points;
	} else {
	  ++cluster_id;
	  std::cerr << "# New cluster: " << cluster_id << std::endl; 
	  (it->second).second = cluster_id;
	  ++visited_points;
	  for (std::size_t n = 0; n < near_points.size(); ++n) {
	    cluster_t::iterator jit = db_map.find( near_points[n]);
	    if (jit != db_map.end() ) {
	      if ((jit->second).second < 0 ) { // not visited
		point_data_t j_near_points;
		kdt.k_nearest_neighbors_in(j_near_points, near_points[n], 1000, knn_maxradius); 
		if (j_near_points.size() >= min_points ) {
		  near_points.insert( near_points.end(), j_near_points.begin(), j_near_points.end());
		}
		(jit->second).second = cluster_id;
		++visited_points;
	      }
	    }
	  }
	}
      }
    }
    std::cerr<< "### DBSCAN: visited points = "<< visited_points << std::endl;
    std::cerr << "### DBSCAN: Finished clustering with " << cluster_id <<" clusters" << std::endl;

    // Write in the cluster output vector
    for (cluster_t::iterator it = db_map.begin(); it != db_map.end(); ++it) {
      id_[it->second.first] = it->second.second;
    }

    fill_cluster_barycenters(cluster_id+1, false);
    compute_cluster_info(cluster_info_, cluster_centers_.size(), false );    
  }

  float cluster_analysis::kmeans_assign(const cluster_analysis::point_data_t& centers) {

    std::size_t data_count = data_.size();
    std::size_t changed_assign = 0;
    std::size_t K = centers.size();
    for (std::size_t s = 0; s < data_count; ++s) {
      point_t csi = data_[s].first;
      int k_assign = -1;
      float min_dist = 1.0e10;
      if (csi[0] != -1.0f ) {
	for (std::size_t k = 0; k < K; ++k) {
	  float d  = csi.distance_to( centers[k]);
	  if (d < min_dist ) {
	    min_dist = d;
	    k_assign = k;
	  }
	}
      }
      if (k_assign != id_[s] ) {
	changed_assign++;
	id_[s] = k_assign;
      }
    }
    return float(changed_assign)/float(data_count);
  }

  void cluster_analysis::kmeans_centroid(cluster_analysis::point_data_t& centers) {
    
    typedef sl::vector3f vector_t;

    std::size_t K = centers.size();
    std::vector< std::pair<vector_t,int> > acc(K);
    for (std::size_t k =0; k < K; ++k) {
      acc[k] = std::make_pair(vector_t(), 0);
    }

    // Accumulate parameters
    std::size_t data_count = data_.size();
    for(std::size_t s = 0; s < data_count; ++s) {
      point_t csi = data_[s].first;
      if (csi[0] != -1.0f) {
	int k = id_[s];
	acc[k].first += as_vector( csi);
	acc[k].second++;
      }
    }
    // Reassign centers
    for (std::size_t k=0; k<K;++k) {
      vector_t c = acc[k].first;
      if (acc[k].second != 0 ) {
	c *= (1.0f/acc[k].second);
      }
      centers[k] = as_point(c);
    }
  }

  void cluster_analysis::kmeans(std::size_t K, cluster_analysis::point_data_t& centers) {

    sl::random::uniform<float> rng;
    
    // If not initialized random choice 
    std::size_t data_count = data_.size();
    id_.resize(data_count);
    if (centers.size() < K ) {
      for (std::size_t s = 0; s < data_count; ++s) {
	id_[s]  = int( rng.value() * K );
      }
      centers.resize(K);
      kmeans_centroid(centers);
    } 

    float change_percent = kmeans_assign(centers); 
    std::cerr <<"Initial change percent == "<< change_percent <<std::endl;
    const float change_threshold = 0.001f;
    std::size_t max_iteration_count = 10000;

    std::size_t it;
    for (it = 0; (change_percent > change_threshold) && it < max_iteration_count; ++it) {
      kmeans_centroid(centers);
      change_percent = kmeans_assign(centers);  
    }
    std::cerr <<"Final change percent == "<< change_percent <<std::endl;
    std::cerr <<"## Converged after "<< (it+1) <<" iterations"<<std::endl;
    fill_cluster_barycenters(K);	
  }
  

  void cluster_analysis::fill_cluster_barycenters(std::size_t k, bool knn) {

    typedef sl::vector3f vector_t;

    std::vector<vector_t> centers(k,vector_t());
    std::vector<std::size_t> count(k,0);
    
    std::size_t data_count = data_.size();
    std::size_t first_valid_id = knn?0:1;
    for (std::size_t s = 0; s < data_count; ++s) {
      int cur_id = id_[s];
      if ( cur_id >= int(first_valid_id) ) {
	count[cur_id]++;
	centers[cur_id] +=  as_vector(data_[cur_id].first);
      }
    }

    cluster_centers_.clear();
    cluster_centers_.resize(k);
    for (std::size_t i = first_valid_id; i < k; ++i) {
      point_t c_i = as_point( 1.0f / sl::max(float(count[i]),1.0f) * centers[i]);
      cluster_centers_[i] = c_i; 
    }
    std::cerr << "#### Computed cluster barycenters " << std::endl;  
  }

  void cluster_analysis::compute_cluster_info(std::vector<cluster_info>& ci, std::size_t k, bool knn) {
    
    ci.clear();
    ci.resize(k);   
    
    // Analyze params
    std::size_t data_count = data_.size();
    int threshold = knn?-1:0;
    float tot_data = 0.0f;
    for (std::size_t s = 0; s < data_count; ++s) {
      int cur_id = id_[s];  
      if ( cur_id > threshold ) {
	tot_data += 1.0f;
	point_t p = data_[s].first;
	value_t v = data_[s].second;
	ci[cur_id].count_ += 1.0f;
	ci[cur_id].pmean_ += as_vector(p);
	ci[cur_id].vmean_ += v;
	for(std::size_t i = 0 ; i < 3; ++i ) {
	  if ( ci[cur_id].pmin_[i] > p[i] )  ci[cur_id].pmin_[i] = p[i];
	  if ( ci[cur_id].pmax_[i] < p[i] )  ci[cur_id].pmax_[i] = p[i];
	}
	if (ci[cur_id].vmin_ > v )  ci[cur_id].vmin_ = v;
	if (ci[cur_id].vmax_ < v )  ci[cur_id].vmax_ = v;
      }
    }

    for (std::size_t c = (knn?0:1); c < k; ++c) {
      value_t w = 1.0f/ci[c].count_;
      ci[c].pmean_ = point_t(w*ci[c].pmean_[0],w*ci[c].pmean_[1],w*ci[c].pmean_[2]);
      ci[c].vmean_ *= w;
      ci[c].count_ /= tot_data;
      ci[c].volume_ = (ci[c].pmax_[1] -  ci[c].pmin_[1])*(ci[c].pmax_[2] -  ci[c].pmin_[2]) 
	* (ci[c].pmax_[0] -  ci[c].pmin_[0]);
    }

    std::cerr<<" ##################### " << std::endl;
    for (std::size_t c = (knn?0:1); c < k; ++c) {
      std::cerr <<" ###  CLUSTER [" << c <<"]: count "<< ci[c].count_*tot_data <<  std::endl;
      std::cerr <<" ###  pmin = ( "<< ci[c].pmin_[0]<<" , "<< ci[c].pmin_[1] << " , "<< ci[c].pmin_[2] << ")" <<  std::endl;
      std::cerr <<" ###  pmax = ( "<< ci[c].pmax_[0]<<" , "<< ci[c].pmax_[1] << " , "<< ci[c].pmax_[2] << ")" <<  std::endl;
      std::cerr <<" ###  pmean = ( "<< ci[c].pmean_[0]<<" , "<< ci[c].pmean_[1] << " , "<< ci[c].pmean_[2] << ")" <<  std::endl;
      float volume = ci[c].volume_;
      std::cerr<<" ### volume: " << volume << ", density: "<< ci[c].count_ /volume << std::endl;
      std::cerr<<" ### vmin: " << ci[c].vmin_  << " , vmax: "<< ci[c].vmax_ << " , vmean: " << ci[c].vmean_ << std::endl;
    }

    std::cerr<<" ##################### " << std::endl;
    std::ofstream out_file("cluster_info.dat");
    out_file <<"ID COUNT VOLUME CX CY CZ VMIN VMAX VMEAN "<< std::endl;
    for (std::size_t c = (knn?0:1); c < k; ++c) {
      out_file<<c<<" "<<ci[c].count_*tot_data<<" "<< ci[c].volume_<<" "
	     <<ci[c].pmean_[0]<<" "<<ci[c].pmean_[1]<<" "<<ci[c].pmean_[2]<<" "
	     <<ci[c].vmin_<<" "<<ci[c].vmax_<<" "<<ci[c].vmean_<<std::endl;
    }    

    out_file.close();
  }


}
