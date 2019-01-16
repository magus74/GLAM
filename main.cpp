//#include <getopt.h>             /* getopt_long() */
#include <string>
#include <iostream>
#include <influence_map_builder.hpp>
#include <parser.hpp>
#include <sl/clock.hpp>
#include <cluster_analysis.hpp>
#include <vic/ply_io/std_triangle_mesh.hpp>
#include <vic/ply_io/std_ply_io.hpp>

#include <boost/program_options.hpp>

static const double M_PI = 3.1415;

namespace po = boost::program_options;

template <class type>
static bool parse_argument(type& val, const po::variables_map& vm,
	const std::string& key) {

	bool result;
	if (vm.count(key)) {
		val = vm[key].as<type>();
		std::cerr << key.c_str() << " set to "
			<< val << std::endl;
		result = true;
	}
	else {
		std::cerr << "MSG: " << key.c_str() << " not set." << std::endl;
		result = false;
	}
	return result;
}

#if 0
static const char short_options [] = "o:s:b:w:r:l:c:t:p:d:f:x:m:e:h";

static const struct option
long_options [] = {
        { "output-file",             required_argument,      0,           'o' },
        { "source-file",                required_argument,      0,           's' },
	{ "laplacian",               no_argument,              0,           'b'},
	{ "weight-max",              required_argument,            0,           'w' },
	{ "influence-radius",        required_argument,            0,           'r' },
	{ "light-model",        required_argument,            0,           'l' },
	{ "colormap",        required_argument,            0,           'c' },
	{ "threshold-cluster",        required_argument,            0,           't' },
	{ "peak-file", required_argument,            0,           'p' },
	{ "cluster-datafile", required_argument,            0,           'd' },
	{ "object-datafile", required_argument,            0,           'f' },
	{ "colorize-per-object", required_argument,            0,           'x' },
	{ "maximum-expected-absorption", required_argument,            0,           'e' },
	{ "help",                    no_argument,            0,           'h' },
        { 0, 0, 0, 0 }
};

static void usage(int /*argc*/, char *argv[]) {
  std::cerr <<
    "Usage: " << argv[0] << " [options] in_ply_filename\n\n" 
    "Options:\n"
    "-o | --output-file <out_ply_filename>          Output file name\n"
    "-s | --source-file <src_filename>              Source file name\n"
    "-b | --laplacian (Beltrami)                    Compute laplacian field\n"
    "-w | --weight-max                              Set weight max \n"
    "-r | --influence-radius                        Set influence radius \n"
    "-l | --light-model                             Set light model (LAMBERT or GAUSSIAN) \n"
    "-c | --colormap                                Set colormap ascii file (t_norm R_byte G_byte B_byte ) \n"
    "-t | --cluster-threshold                       Set normalized threshold for clustering \n"
    "-p | --peak-file                               Output cluster peak ply filename \n"
    "-d | --cluster-datafile                        Output data filename \n"
    "-f | --object-datafile                         Output object data filename \n"
    "-x | --colorize-per-object                     Output cluster peak ply filename \n"
    "-m | --colormap-per-object                     Set colormap ascii file (t_norm R_byte G_byte B_byte ) \n"
    "-e | --maximum-expected-absorption             Set maximum expected absorption \n"
    "-h | --help                                    Print this message\n"
    "";
}

#endif


typedef std::vector< std::pair<float,vic::influence_map_builder::color_t> > colormap_t;

bool import_colormap(colormap_t& cmap, const std::string& colormap_file) {

  bool result = true;
  std::ifstream fin( colormap_file.c_str());
  cmap.clear();
  if ( fin.is_open() ) {
    while ( !fin.eof() ) {
      float t;
      fin >> t;
      int r;
      fin >> r;
      int g;
      fin >> g;
      int b;
      fin >> b;
      std::cerr << t << " "<< r <<" " << g <<" "<<b<<" "<< std::endl;
      cmap.push_back( std::make_pair(t, vic::influence_map_builder::color_t(float(r)/255.0f,float(g)/255.0f,float(b)/255.0f)));
    }
  } else {
    result = false;
    std::cerr<<"### Error importing color map from file " << colormap_file << std::endl;
  }
  return result;
}
 
// Jet colormap from blue to red 
void fill_default_colormap(colormap_t& cmap) {
  cmap.clear();
#if 1
  cmap.push_back( std::make_pair(0.0f, vic::influence_map_builder::color_t(0.0f,0.0f,1.0f)));
  cmap.push_back( std::make_pair(0.25f, vic::influence_map_builder::color_t(0.0f,1.0f,1.0f)));
  cmap.push_back( std::make_pair(0.5f, vic::influence_map_builder::color_t(0.0f,1.0f,0.0f)));
  cmap.push_back( std::make_pair(0.75f, vic::influence_map_builder::color_t(1.0f,1.0f,0.0f)));
  cmap.push_back( std::make_pair(1.0f, vic::influence_map_builder::color_t(1.0f,0.0f,0.0f)));
#else
  cmap.push_back( std::make_pair(0.0f, vic::influence_map_builder::color_t(237.0f/255.0f,248.0f/255.0f,251.0f/255.0f)));
  cmap.push_back( std::make_pair(0.25f, vic::influence_map_builder::color_t(179.0f/255.0f,205.0f/255.0f,227.0f/255.0f)));
  cmap.push_back( std::make_pair(0.5f, vic::influence_map_builder::color_t(140.0f/255.0f,150.0f/255.0f,198.0f/255.0f)));
  cmap.push_back( std::make_pair(0.75f, vic::influence_map_builder::color_t(136.0f/255.0f,86.0f/255.0f,167.0f/255.0f)));
  cmap.push_back( std::make_pair(1.0f, vic::influence_map_builder::color_t(129.0f/255.0f,15.0f/255.0f,124.0f/255.0f)));
#endif
}


void push_sphere_to(vic::std_triangle_mesh& mesh,
		    const sl::point3f& gamma,
		    const vic::influence_map_builder::color_t& col,
		    float rad,
		    float q = 0.0f,
		    std::size_t sampling = 16) {
  
  std::size_t  start_index = mesh.vertex_count();
  //  std::cerr << "start index " << start_index << std::endl; 

  for (std::size_t i = 0; i < sampling-1; ++i ) {
    for (std::size_t j = 0; j < sampling-1; ++j ) {
      sl::tuple2i i00(i,j);
      sl::point2f a00(-M_PI + float(i00[0])/float(sampling-1) * 2.0f * M_PI,
			    -0.5f*M_PI + float(i00[1])/float(sampling-1)  * M_PI);
      sl::point3f p00 = gamma + rad * sl::vector3f( cosf(a00[0])*cosf(a00[1]), sinf(a00[0])*cosf(a00[1]), sinf(a00[1]));

      sl::tuple2i i10(i+1,j);
      sl::point2f a10(-M_PI + float(i10[0])/float(sampling-1) * 2.0f * M_PI,
		      -0.5f*M_PI + float(i10[1])/float(sampling-1)  * M_PI);
      sl::point3f p10 = gamma + rad * sl::vector3f( cosf(a10[0])*cosf(a10[1]), sinf(a10[0])*cosf(a10[1]), sinf(a10[1]));

      sl::tuple2i i11(i+1,j+1);
      sl::point2f a11(-M_PI + float(i11[0])/float(sampling-1) * 2.0f * M_PI,
			    -0.5f*M_PI + float(i11[1])/float(sampling-1)  * M_PI);
      sl::point3f p11 = gamma + rad * sl::vector3f( cosf(a11[0])*cosf(a11[1]), sinf(a11[0])*cosf(a11[1]), sinf(a11[1]));

      sl::tuple2i i01(i,j+1);
      sl::point2f a01(-M_PI + float(i01[0])/float(sampling-1) * 2.0f * M_PI,
			    -0.5f*M_PI + float(i01[1])/float(sampling-1)  * M_PI);
      sl::point3f p01 = gamma + rad * sl::vector3f( cosf(a01[0])*cosf(a01[1]), sinf(a01[0])*cosf(a01[1]), sinf(a01[1]));

      std::size_t id00 = start_index + i     + j*sampling;
      std::size_t id10 = start_index + i + 1 + j*sampling;
      std::size_t id11 = start_index + i + 1 + (j+1)*sampling;
      std::size_t id01 = start_index + i     + (j+1)*sampling;

      mesh.push_vertex( p00, sl::point3f(col[0],col[1],col[2]),
			vic::std_triangle_mesh::dual_vector_t(cosf(a00[0])*cosf(a00[1]), sinf(a00[0])*cosf(a00[1]), sinf(a00[1])),q);
      mesh.push_vertex( p10, sl::point3f(col[0],col[1],col[2]),
			vic::std_triangle_mesh::dual_vector_t(cosf(a10[0])*cosf(a10[1]), sinf(a10[0])*cosf(a10[1]), sinf(a10[1])),q);
      mesh.push_vertex( p11, sl::point3f(col[0],col[1],col[2]),
			vic::std_triangle_mesh::dual_vector_t(cosf(a11[0])*cosf(a11[1]), sinf(a11[0])*cosf(a11[1]), sinf(a11[1])),q);
      mesh.push_vertex( p01, sl::point3f(col[0],col[1],col[2]),
			vic::std_triangle_mesh::dual_vector_t(cosf(a01[0])*cosf(a01[1]), sinf(a01[0])*cosf(a01[1]), sinf(a01[1])),q);

      std::size_t idx = mesh.vertex_count();

      mesh.push_triangle( vic::std_triangle_mesh::tri_t(idx-4,idx-3,idx-2));
      mesh.push_triangle( vic::std_triangle_mesh::tri_t(idx-2,idx-1,idx-4));

    }
  }
  
  //  std::cerr << "end index " << mesh.vertex_count() << std::endl; 
}



typedef std::vector<vic::cluster_info> cluster_info_vector_t;

static void export_peaks_info(const std::vector< cluster_info_vector_t>& ci, 
			      const std::string& out_file) {

  std::ofstream fout(out_file.c_str());
  fout<<"Obj VertNormCount Volume Vmax Vmin Vavg" << std::endl;
  for (std::size_t i = 1; i < ci.size(); ++i ) {
    std::vector<vic::cluster_info> cl_info = ci[i];
    for (std::size_t j = 1; j <  cl_info.size(); ++j ) {
      fout<<i<<" "<<cl_info[j].count_<<" "<<cl_info[j].volume_ << " "<<cl_info[j].vmax_<<" "<<cl_info[j].vmin_<<" "<<cl_info[j].vmean_<<std::endl;
    }
  }
  fout.close();
  std::cerr<<"### Exported cluster info in file "<< out_file << std::endl;
}

static void export_object_info(const std::vector< vic::object_info> & obj_info,
			       const std::string& out_file) {

 std::ofstream fout(out_file.c_str());
 fout<<"Obj Amax Atot Stot Fmax Ftot xmin ymin zmin xmax ymax zmax" << std::endl;
 for (std::size_t i =0; i < obj_info.size(); ++i ) {
   vic::object_info info = obj_info[i];
   fout<<info.id_<<" "<<info.alpha_max_<<" "<<info.alpha_tot_
       <<" "<<info.sigma_tot_<<" "<<info.phi_max_<<" "<<info.phi_tot_
       <<" "<<info.box_[0][0]<<" "<<info.box_[0][1]<<" "<<info.box_[0][2]
       <<" "<<info.box_[1][0]<<" "<<info.box_[1][1]<<" "<<info.box_[1][2] << std::endl;
 }
 fout.close();
  std::cerr<<"### Exported object info in file "<< out_file << std::endl;
}

static void export_obj(const vic::std_triangle_mesh& mesh,
		       const std::string& out_file) {

  std::ofstream of(out_file.c_str());
  of <<"mtllib peak.mtl"<<std::endl;
  of <<"usemtl mat"<<std::endl;
  
  for(std::size_t i = 0; i < mesh.vertex_count(); ++i) {
    vic::std_triangle_mesh::point_t p = mesh.get_position(i);
    vic::std_triangle_mesh::dual_vector_t n = mesh.get_normal(i);
    vic::std_triangle_mesh::value_t q =  mesh.get_quality(i);
    of << "v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
    of << "vt "<<q<<" "<<"0.5 "<<std::endl;
    of << "vn "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
  }
  for(std::size_t i = 0; i < mesh.triangle_count(); ++i ) {
    vic::std_triangle_mesh::tri_t tr = mesh.get_triangle(i);
    of <<"f "<<(tr[0]+1)<<"/"<<(tr[0]+1)<<"/"<<(tr[0]+1)
       <<" "<<(tr[1]+1)<<"/"<<(tr[1]+1)<<"/"<<(tr[1]+1)
       <<" "<<(tr[2]+1)<<"/"<<(tr[2]+1)<<"/"<<(tr[2]+1)<< std::endl;
  }
  of.close();
}


static void export_cluster_peaks(const std::vector<cluster_info_vector_t >& cluster_info, 
				 const colormap_t& cmap,
				 const std::string& out_file,
				 float vmin,
				 float wmax) {

  vic::std_triangle_mesh mesh;
  sl::interpolation_track<float,vic::influence_map_builder::color_t> color_map;
  color_map.set_continuity(0);
  for(std::size_t i=0; i< cmap.size(); ++i) {
    color_map.insert( cmap[i]);
  }     

  cluster_info_vector_t ci;
  for (std::size_t i=0; i< cluster_info.size(); ++i) {
    cluster_info_vector_t::const_iterator cluster_info_first = ++cluster_info[i].begin();
    ci.insert(ci.end(),cluster_info_first,cluster_info[i].end());
  }

  std::vector<float> importance(ci.size()); 
  float importance_max = 0.0f;
  float vmax = 0.0f;
  float rmax = 0.0f;
  for (std::size_t i = 0; i < ci.size(); ++i ) {
    importance[i] = ci[i].volume_ * ci[i].vmean_;
    rmax = sl::max(rmax, ci[i].pmin_.distance_to( ci[i].pmax_));
    if ( importance_max < importance[i]) importance_max = importance[i];
    if ( vmax < ci[i].vmax_ ) vmax = ci[i].vmax_;
  }
  
  if (wmax != 0 ) {
    vmax = wmax; 
  }

  // Draw a sphere for each cluster
  for (std::size_t i = 0; i < ci.size(); ++i ) {
    float q = (ci[i].vmax_-vmin)/(vmax-vmin);
#if 1    
    vic::influence_map_builder::color_t col = color_map.value_at (q);
    float r_i = 0.25f * ci[i].pmin_.distance_to( ci[i].pmax_);
    std::cout << (r_i*ci[i].vmean_/vmax) << std::endl;
    float rad = sl::median(0.02f, r_i*ci[i].vmean_/vmax, 0.2f);
#else
    vic::influence_map_builder::color_t col = color_map.value_at ( (ci[i].vmax_ - vmin)/(vmax - vmin));
    float rad = 0.5f*( ci[i].pmin_.distance_to( ci[i].pmax_) ); //ci[i].vmax_ / vmax * rmax;
#endif
    //std::cerr <<" Cluster radius "<< rad << std::endl;
    push_sphere_to( mesh, ci[i].pmean_, col, rad, q);
  }
  if ( out_file.find(".ply") != std::string::npos ) {
    vic::std_ply_io::to_file(mesh, out_file.c_str());
  } else if ( out_file.find(".obj") != std::string::npos ) {
    export_obj(mesh, out_file);
  } else {
    std::cerr<<" ##### error: wrong output format in file " << out_file << std::endl; 
    return;
  }
  std::cerr <<"# Written file " << out_file << std::endl;
}


int main(int argc, char** argv) {


	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("input,i", po::value<std::string>(), "set input file (ply) ")
		("source,s", po::value<std::string>(), "set source file (ply)")
		("laplacian,b", po::bool_switch(), "compute laplacian field (default false)")
		("weight,w", po::value<float>(), "set weight max")
		("radius,r", po::value<float>(), "set influence radius (default 0.5)")
		("lambertian,l", po::bool_switch(), "set lambertian light model (default false = gaussian)")
		("colormap,c", po::value<std::string>(), "set colormap ascii file (t_norm R_byte G_byte B_byte )")
		("threshold,t", po::value<float>(), "set normalized threshold for clustering")
		("peaks,p", po::value<std::string>(), "output cluster peak ply file")
		("cluster-data,d", po::value<std::string>(), "output cluster data file")
		("object-data,f", po::value<std::string>(), "output object data filename")
		("colorized,x", po::value<std::string>(), "output cluster peak ply filename ")
		("colormap-per-object,m", po::value<std::string>(), "colormap per object ascii file (t_norm R_byte G_byte B_byte ) ")
		("expected-absorption,e", po::value<float>(), "set maximum expected absorption")
		("output,o", po::value<std::string>(), "output file ply")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	po::notify(vm);

	if (vm.count("help")) {
		std::cerr << desc << "\n";
		return 1;
	}

  std::string in_ply_fname;
  std::string out_ply_fname;
  std::string src_fname;
  std::string colormap_file;
  std::string object_colormap_file;
  bool lambert_light_model = false;
  bool laplacian_field = false;
  float wmax = 0.0f;
  float radius = 0.5f;
  float cluster_threshold = 0.0f;
  float amax = 0.0f; 
  std::string peak_file = "glam_peaks.ply";
  std::string cluster_data_file = "cluster_data.dat";
  std::string object_data_file = "object_data.dat";
  std::string object_ply_file = "objects.ply";
 
  // Parse mandatory arguments
  bool valid_arg = true;
  valid_arg = parse_argument(in_ply_fname, vm, "input");
  if (!valid_arg) {
	  std::cerr << " Exiting " << std::endl;
	  return 1;
  }

  valid_arg = parse_argument(out_ply_fname, vm, "output");
  if (!valid_arg) {
	  std::cerr << " Exiting " << std::endl;
	  return 1;
  }

  valid_arg = parse_argument(src_fname, vm, "source");
  if (!valid_arg) {
	  std::cerr << " Exiting " << std::endl;
	  return 1;
  }

  // Parse optional arguments
  valid_arg = parse_argument(colormap_file, vm, "colormap");
  valid_arg = parse_argument(object_colormap_file, vm, "colormap-per-object");
  valid_arg = parse_argument(lambert_light_model, vm, "lambertian");
  valid_arg = parse_argument(laplacian_field, vm, "laplacian");
  valid_arg = parse_argument(wmax, vm, "weight");
  valid_arg = parse_argument(radius, vm, "radius");
  valid_arg = parse_argument(cluster_threshold, vm, "threshold");
  valid_arg = parse_argument(amax, vm, "expected-absorption");
  valid_arg = parse_argument(peak_file, vm, "peaks");
  valid_arg = parse_argument(cluster_data_file, vm, "cluster-data");
  valid_arg = parse_argument(object_data_file, vm, "object-data");
  valid_arg = parse_argument(object_ply_file, vm, "colorized");


#if 0
  for (;;) {
    int index;
    int c;
                
    c = getopt_long(argc, argv,
		    short_options, long_options,
		    &index);

    if (-1 == c)
      break;
    
    switch (c) {
    case 0: /* getopt_long() flag */
      break;
      
    case 'o':
      out_ply_fname=std::string(optarg);
      break;
    case 's':
      src_fname=std::string(optarg);
      break;
    case 'b':
      laplacian_field = true;
      break;
    case 'w':
      wmax = float(std::atof(optarg));
      break;
    case 'r':
      radius = float(std::atof(optarg));
      break;
    case 'l':
      light_model = std::string(optarg);
      break;
    case 'c':
      colormap_file = std::string(optarg);
      break;
    case 't':
      cluster_threshold = float(std::atof(optarg));
      break;
    case 'p':
      peak_file = std::string(optarg);  
      break;
    case 'd':
      cluster_data_file = std::string(optarg);
      break;
    case 'f':
      object_data_file = std::string(optarg);
      break;
    case 'x':
      object_ply_file = std::string(optarg);
      break;
    case 'm':
      object_colormap_file = std::string(optarg);
      break;
    case 'e':
      amax = float(std::atof(optarg));
      break;
    case 'h':
      usage (argc, argv);
      exit (0);
      break;
    default:
      usage (argc, argv);
      exit (1);
      break;
    }
  }
  
  //  for (index = optind; index < argc; index++)
  if (argc != optind+1) {
    std::cerr << "Missing arguments" << std::endl;
    usage (argc, argv);
    exit (1);
  }

  std::string in_ply_fname=argv[optind];
  if (out_ply_fname.empty()) {
    out_ply_fname=sl::pathname_without_extension(in_ply_fname) + "_out.ply";
  }

  if (src_fname.empty()) {
    std::cerr << "Missing source filename" << std::endl;
    exit (1);
  }
#endif

  vic::parser src_parser(src_fname); 
  vic::influence_map_builder builder(src_parser);
  builder.load_model(in_ply_fname);

  if (wmax > 0.0f ) {
    builder.set_max_weight(wmax);
  } else {
    builder.set_update_max_weight(true);
  }
  if (radius > 0.0f ) {
    builder.set_max_influence_distance(radius);
  }

  if (lambert_light_model ) {
    builder.set_light_model( vic::influence_map_builder::LAMBERT );
  }

  std::vector< std::pair<float,vic::influence_map_builder::color_t> > cmap;
  bool ok_import_colormap = false;
  if ( colormap_file.size() > 0 ) {
    ok_import_colormap = import_colormap(cmap, colormap_file);
  }
  if ( ok_import_colormap ) {
    std::cerr << "### Imported color map from file " << colormap_file << std::endl;
  } else {
    fill_default_colormap(cmap);
  }

  sl::real_time_clock ck;
  builder.compute_color_map_field( cmap, laplacian_field, true);
  float proc_time_sec = ck.elapsed().as_milliseconds() / 1000.0f; 
  std::cerr<<"#### Proc T  == " << proc_time_sec << " s " << std::endl;
  double surface_flux = builder.compute_surface_flux();
  std::cerr<<"  ## Expected absorption rate  == " <<  surface_flux  << std::endl;
  export_object_info( builder.obj_info(), object_data_file );


  if (cluster_threshold > 0.0f ) {
    std::vector < vic::influence_map_builder::peak_vector_t > peaks;
    builder.export_peak_positions( cluster_threshold, peaks);
    std::vector< cluster_info_vector_t> cluster_data;
    for ( std::size_t i = 0; i < peaks.size(); ++i ) {
      vic::cluster_analysis cl(peaks[i]);
      cl.dbscan(0.2f*radius,2);
      cluster_data.push_back( cl.cluster_info_data());
    }
    colormap_t cmap;
    fill_default_colormap(cmap);
    export_cluster_peaks( cluster_data, cmap, peak_file, cluster_threshold*wmax,  wmax );
    export_peaks_info( cluster_data, cluster_data_file);
  }

  builder.save_model( out_ply_fname);

  if ( object_colormap_file.size() > 0 ) {
    import_colormap(cmap, object_colormap_file);
  }
  builder.colorize_objects( cmap, amax, (cluster_threshold > 0.0f? cluster_threshold : 0.1f));
  builder.save_model( object_ply_file);
  

  std::cout<<"### Finished" << std::endl;
  return 0;
}
