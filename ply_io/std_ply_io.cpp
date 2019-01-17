#include <vic/ply/ply.h>
#include <vic/ply_io/std_ply_io.hpp>
#include <stdio.h>

#define MAT2IDX(i,j,col) i * col + j

namespace vic {

  /* user's vertex and face definitions for a polygonal object */
  typedef vic::std_triangle_mesh::uint_t                            uint_t;
  typedef sl::point3f                                               point_t;
  typedef sl::point3f                                               color_t;
  typedef sl::row_vector3f                                          dual_vector_t;
  typedef sl::fixed_size_point<3,uint_t>                            tri_t;
  typedef unsigned char uchar;
  typedef sl::fixed_size_vector<sl::column_orientation,4,float>     hposition_t;      

  typedef struct Vertex {
    float x,y,z;
    float nx,ny,nz;
    unsigned char red,green,blue;
    float quality;
    void *other_props;       /* other properties */
  } Vertex;
      
  typedef struct Face {
    unsigned char nverts;    /* number of vertex indices in list */
    int *verts;              /* vertex index list */
    void *other_props;       /* other properties */
  } Face;

  static PlyProperty vert_props[] = { /* list of property information for a vertex */
    {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
    {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
    {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
    {"nx", Float32, Float32, offsetof(Vertex,nx), 0, 0, 0, 0},
    {"ny", Float32, Float32, offsetof(Vertex,ny), 0, 0, 0, 0},
    {"nz", Float32, Float32, offsetof(Vertex,nz), 0, 0, 0, 0},
    {"red",Uint8,Uint8,offsetof(Vertex,red), 0, 0, 0, 0},
    {"green",Uint8,Uint8,offsetof(Vertex,green), 0, 0, 0, 0},
    {"blue",Uint8,Uint8,offsetof(Vertex,blue), 0, 0, 0, 0},
    {"quality",Float32,Float32,offsetof(Vertex,quality), 0, 0, 0, 0},
  };

  static PlyProperty face_props[] = { /* list of property information for a face */
    {"vertex_indices", Int32, Int32, offsetof(Face,verts),
     1, Uint8, Uint8, offsetof(Face,nverts)},
  };


  void std_ply_io::to_file(const std_triangle_mesh& tm, const char* str) {
    //-------------------------------//
    PlyFile* out_ply = 0;
    int nelems = 2;
    uint_t nverts,nfaces;
    char *elem_names[] = {"vertex","face"};
    out_ply = open_for_writing_ply((char *)str,nelems,elem_names,PLY_BINARY_LE);

    nverts = tm.vertex_count();
    nfaces = tm.triangle_count();

    element_count_ply (out_ply, "vertex", nverts);
    ply_describe_property (out_ply, "vertex", &vert_props[0]);
    ply_describe_property (out_ply, "vertex", &vert_props[1]);
    ply_describe_property (out_ply, "vertex", &vert_props[2]);
    ply_describe_property (out_ply, "vertex", &vert_props[3]);
    ply_describe_property (out_ply, "vertex", &vert_props[4]);
    ply_describe_property (out_ply, "vertex", &vert_props[5]);
    ply_describe_property (out_ply, "vertex", &vert_props[6]);
    ply_describe_property (out_ply, "vertex", &vert_props[7]);
    ply_describe_property (out_ply, "vertex", &vert_props[8]);
    ply_describe_property (out_ply, "vertex", &vert_props[9]);
	
    element_count_ply (out_ply, "face", nfaces);
    ply_describe_property (out_ply, "face", &face_props[0]);

    header_complete_ply (out_ply);

    put_element_setup_ply (out_ply, "vertex");
    Vertex* v_ptr = new Vertex;
    for(uint_t i = 0; i < nverts; ++i) {
      if (i%100000 == 0) std::cout << "Writing vertex " << i+1 << "/" << nverts << " - (" << (i+1)*100/nverts << "%)" << "\r";
      const point_t p = tm.get_position(i);
      const color_t c = tm.get_color(i);
      const dual_vector_t n = tm.get_normal(i);
      const float q  = tm.get_quality(i);
      v_ptr->x = p[0];
      v_ptr->y = p[1];
      v_ptr->z = p[2];
      v_ptr->nx = n[0];
      v_ptr->ny = n[1];
      v_ptr->nz = n[2];
      v_ptr->red = uchar(c[0]*255);
      v_ptr->green = uchar(c[1]*255);
      v_ptr->blue = uchar(c[2]*255);
      v_ptr->quality = q;
      put_element_ply (out_ply, (void *) v_ptr);
    }
    std::cout << std::endl;
    std::cout << std::endl;
    put_element_setup_ply (out_ply, "face");
    Face* f_ptr = new Face;
    f_ptr->verts = new int[3];
    for(uint_t i = 0; i < nfaces; ++i) {
      if (i%100000 == 0) std::cout << "Writing faces " << i+1 << "/" << nfaces << " - (" << (i+1)*100/nfaces << "%)" << "\r";
      const tri_t t = tm.get_triangle(i);
      f_ptr->nverts = 3;
      f_ptr->verts[0] = t[0];
      f_ptr->verts[1] = t[1];
      f_ptr->verts[2] = t[2];
      put_element_ply (out_ply, (void *) f_ptr);
    }
    std::cout << std::endl;
    std::cout << std::endl;
    fclose(out_ply->fp);
    free(out_ply);
  }

  void std_ply_io::from_file(const char* str, std_triangle_mesh& tm, float scale) {
    //    std::cout << "____TEST____" << std::endl;
    PlyFile* in_ply = 0;
    FILE* f = fopen(str, "r");  
    if (f == NULL) {
      std::cout << "problems occur opening file" << std::endl;
    } else {
      std::cout << "file opened" << std::endl;
    }
    if (f) { 
      in_ply = read_ply(f);
    }
    	      
    if (in_ply) {
      bool corrupted_file  = false;
	
      uint_t nverts,nfaces;
	
      uint_t i,j;
      int elem_count;
      char *elem_name;
	
      /* examine each element type that is in the file (vertex, face) */
      for (i = 0; i < (uint_t)in_ply->num_elem_types && !corrupted_file; i++) {
	/* prepare to read the i'th list of elements */
	elem_name = setup_element_read_ply (in_ply, i, &elem_count);
	  
	if (equal_strings ("vertex", elem_name)) {
	  nverts = elem_count;
	  tm.resize_vertex(nverts);
	    
	  /* set up for getting vertex elements */
	  /* (we want x,y,z) */
	  setup_property_ply (in_ply, &vert_props[0]);
	  setup_property_ply (in_ply, &vert_props[1]);
	  setup_property_ply (in_ply, &vert_props[2]);
	  setup_property_ply (in_ply, &vert_props[3]);
	  setup_property_ply (in_ply, &vert_props[4]);
	  setup_property_ply (in_ply, &vert_props[5]);
	  setup_property_ply (in_ply, &vert_props[6]);
	  setup_property_ply (in_ply, &vert_props[7]);
	  setup_property_ply (in_ply, &vert_props[8]);
	  setup_property_ply (in_ply, &vert_props[9]);
	    
	  /* grab the vertex elements and store them in our list */
	  for (j = 0; j < nverts; j++) {
	    if (j%100000 == 0) std::cout << "Reading vertex " << j+1 << "/" << nverts << " - (" << (j+1)*100/nverts << "%)" << "\r";
	    Vertex v;
	    //default v - usefull if no data was found
	    v.x = 0;
	    v.y = 0;
	    v.z = 0;
	    v.nx = 0;
	    v.ny = 0;
	    v.nz = 0;
	    v.red = 195;
	    v.green = 195;
	    v.blue = 195;
	    v.quality = 0;
	    //
	    get_element_ply (in_ply, (void *) &v);
	    const float MAX_COORD_VALUE = 1e20;
	    const float MIN_COORD_VALUE = -MAX_COORD_VALUE;
	    if (v.x<MIN_COORD_VALUE || v.x>MAX_COORD_VALUE || v.y<MIN_COORD_VALUE || v.y>MAX_COORD_VALUE || v.z<MIN_COORD_VALUE || v.z>MAX_COORD_VALUE) {
	      corrupted_file = true;
	      std::cerr << "possible wrong input vertex coordinate " << v.x << " " << v.y << " " << v.z << std::endl;
	    } else {
	      tm.set_position(j,point_t(v.x * scale, v.y * scale ,v.z * scale));
	      tm.set_color(j,color_t(v.red/255.0f,v.green/255.0f,v.blue/255.0f));
	      tm.set_normal(j,dual_vector_t(v.nx,v.ny,v.nz));
	      tm.set_quality(j,float(v.quality));
	    }
	  }
	} else if (equal_strings ("face", elem_name)) {
	  /* create a list to hold all the face elements */
	  nfaces = elem_count;
	  tm.resize_triangle(nfaces);
	    
	  /* set up for getting face elements */
	  /* (all we need are vertex indices) */
	  setup_property_ply (in_ply, &face_props[0]);
	  std::cout << "\n";
	  /* grab all the face elements and place them in our mesh */
	  for (j = 0; j < nfaces && !corrupted_file; j++) {
	    if (j%100000 == 0) std::cout << "Reading face " << j+1 << "/" << nfaces << " - (" << (j+1)*100/nfaces << "%)" << "\r";
	    Face f;
	    f.nverts = 0;
	    f.verts = 0;
	    get_element_ply (in_ply, (void *) &f);
	    if (!f.verts) {
	      std::cerr << "PLY: skipping empty face with " << (int)f.nverts << " vertices." << std::endl;
	    } else {
#if 1
	      if (f.nverts != 3) {
		std::cerr << "PLY: skipping face with " << (int)f.nverts << " vertices." << std::endl;
	      } else {
		tm.set_triangle(j,tri_t(f.verts[0],f.verts[1],f.verts[2]));
	
		if (f.verts[0] < 0 || f.verts[0] >= nverts ||  
		    f.verts[1] < 0 || f.verts[1] >= nverts || 
		    f.verts[2] < 0 || f.verts[2] >= nverts) {
		  std::cerr << "wrong vertex index in face " << f.verts[0] << " " <<  f.verts[1] << " " <<  f.verts[2] << " not in 0, " << nverts << std::endl;
		  corrupted_file = true;
		} else {

		  //FIXME insert normal computation for each triangle
		  const tri_t t = tm.get_triangle(j);
		  const point_t p0 = tm.get_position(t[0]);
		  const point_t p1 = tm.get_position(t[1]);
		  const point_t p2 = tm.get_position(t[2]);

		  dual_vector_t dv = sl::normal(p0,p1,p2);
#if 1
		  float dv_norm = dv.two_norm();
		  if (dv_norm>1e-5) {
		    dv /= dv_norm;
		    for(std::size_t idx = 0; idx<3; ++idx) {
		      dual_vector_t n0 = tm.get_normal(t[idx]);
		      n0 += dv;
		      tm.set_normal(t[idx],n0);
		    }
		  }
#else
		  float e01_len2 = (p1-p0).two_norm_squared();
		  float e12_len2 = (p2-p1).two_norm_squared();
		  float e20_len2 = (p0-p2).two_norm_squared();
		  float w0 = e01_len2*e20_len2;
		  float w1 = e01_len2*e12_len2;
		  float w2 = e12_len2*e20_len2;

		  if (w0) {
		    dual_vector_t n0 = tm.get_normal(t[0]);
		    n0 += dual_vector_t(float(dv[0]/w0),float(dv[1]/w0),float(dv[2]/w0));
		    tm.set_normal(t[0],n0);
		  }
		  if (w1) {
		    dual_vector_t n1 = tm.get_normal(t[1]);
		    n1 += dual_vector_t(float(dv[0]/w1),float(dv[1]/w1),float(dv[2]/w1));
		    tm.set_normal(t[1],n1);
		  }
		  if (w2) {
		    dual_vector_t n2 = tm.get_normal(t[2]);
		    n2 += dual_vector_t(float(dv[0]/w2),float(dv[1]/w2),float(dv[2]/w2));
		    tm.set_normal(t[2],n2);
		  }
#endif
		} 
	      } 



#else
	      if ((f.nverts != 3) &&  (f.nverts != 4)) {
		std::cerr << "PLY: skipping face with " << (int)f.nverts << " vertices." << std::endl;
	      } else if (f.nverts == 3) {
		tm.set_triangle(j,tri_t(f.verts[0],f.verts[1],f.verts[2]));
		
		//FIXME insert normal computation for each triangle
	      } else {
		tm.resize_triangle(tm.triangle_count()+1);
		nfaces = tm.triangle_count();
		tm.set_triangle(j,tri_t(f.verts[0],f.verts[1],f.verts[3]));
		j+=1;
		tm.set_triangle(j,tri_t(f.verts[1],f.verts[2],f.verts[3]));
	      }
#endif
	      free(f.verts);
	    }
	  }

	  if (!corrupted_file) {
	    //normalizing normal
	    std::cerr << std::endl<< "Normalizing normals" << std::endl;
	    for(uint_t i = 0; i < tm.vertex_count(); ++i) {
	      dual_vector_t n = tm.get_normal(i);
	      float n_len = n.two_norm();
	      if (n_len>1e-4) {
		n /= n_len;
	      } else {
		// std::cerr << "WWW: Null normal set to zero" << std::endl;
		n *= 0.0; // FIXME !!!
	      }
	      tm.set_normal(i,n);
	    }
	  }
	  std::cerr << " ...normalized "<< std::endl;
	} else  {
#if 0
	  /* all non-vertex and non-face elements are grabbed here */
	  get_other_element_ply (in_ply);
#endif
	}
      }
	
      /* close the file */
      free_ply(in_ply);
      fclose(f);
      
      std::cerr << " ...after closing ply"<< std::endl;
      /*FIXME free vlist flist??*/
	
      if (corrupted_file) {
	std::cerr << "Corrupted file " << str  << std::endl;
      }
    } else {
      std::cout << "cannot read ply" << std::endl;
    }
  }

  void std_ply_io::apply_transform(const char* str_dest,
				   const char* str_src,
				   rigid_body_map_t& R,
				   rigid_body_map_t& T) {
    
    std_triangle_mesh tm;
    vic::std_ply_io::from_file(str_src, tm);

    std::size_t nverts = tm.vertex_count();
    for (std::size_t i = 0; i < nverts; ++i) {
      point_t p = tm.get_position(i);
      p = T*R*p;
      tm.set_position(i,p);
    }

    vic::std_ply_io::to_file(tm,str_dest);
  }
  
 
  

} // namespace vic


