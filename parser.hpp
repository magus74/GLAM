#ifndef PARSER_HPP
#define PARSER_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/kdtree.hpp>


namespace vic {

  class parser {
  public:
    typedef sl::point3f point_t;
    typedef sl::kdtree<3,float,point_t> kdtree_t;

  public:
    parser(const std::string& srcfile) {
      create_kdtree(srcfile);
    }

    void nearest_neighbors(std::vector<point_t>& knn,
			   const point_t& p,
			   float max_distance) {
      kdtree_.k_nearest_neighbors_in(knn,
				     p,
				     100000, // Maximum number of neighbors
				     max_distance); // Maximum distance
    }

  protected:
    void create_kdtree(const std::string& srcfile);

  protected:
    kdtree_t  kdtree_;
  };

}

#endif
