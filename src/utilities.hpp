#ifndef BIOUTIL
#define BIOUTIL

// -D_GLIBCXX_ASSERTIONS

#include <vector>
#include <unordered_map>
#include <array>
#include <string>
#include <cmath>
#include <numeric>      // std::iota, std::accumulate
#include <algorithm>    // std::sort, std::stable_sort

constexpr auto PI = 3.14159265358979323846;

struct MCWPARMS {
  std::string common_name;
  std::string extension_route;
  std::string confidence;
  std::string clade;
  double a1;      
  double a2;                    
  double ratio_intercept;   
  double ratio_cr;                
  double ratio_dbh; 
  double ratio_phi;
  double dbh_max_fit;
  double dbh_min_fit; 
  double cw_ceiling;
};

struct Point {
    double x, y;
};

struct TREE {
    double dbh;
    double height;
};


// sorting machinery for generated indices to decreasing sorted vector 
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T> &v, bool increasing = false ) 
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  std::stable_sort(idx.begin(), idx.end(),
       [&v, increasing](size_t i1, size_t i2) { return increasing ? v[i1] < v[i2] : v[i1] >= v[i2]; });

  return idx;
}


double compute_ba( const std::vector<double> &dbh, const std::vector<double> &expansion_factor, const bool imperial );

std::vector<double> compute_bal( const std::vector<double> &dbh, const std::vector<double> &expansion_factor, const bool imperial );

std::vector<double> compute_ccfl( const std::vector<double> &dbh, const std::vector<double> &mcw, const std::vector<double> &expansion, const bool imperial );

double compute_cch( const std::vector<int>    &species,
                    const double ht,
                    const std::vector<double> &dbh,
                    const std::vector<double> &height,
                    const std::vector<double> &crown_length,
                    const std::vector<double> &dacb,        
                    const std::vector<double> &lcw,
                    const std::vector<double> &expansion,
                    const std::unordered_map<int,std::vector<double>> &parameters,
                    const bool imperial );

double compute_dominant_dbh( const std::vector<double> &height,
                             const std::vector<double> &dbh,
                             const std::vector<double> &expansion,
                             const int dominant_cohort_size,
                             const int method );  

double compute_dominant_height( const std::vector<double> &height,
                                const std::vector<double> &dbh,
                                const std::vector<double> &expansion,
                                const int dominant_cohort_size,
                                const int method );                                  

double compute_qmd( const std::vector<double> &dbh,
                    const std::vector<double> &expansion );

double compute_relative_spacing( const std::vector<double> &expansion,
                                 const double dominant_height,
                                 const bool imperial );  

double compute_curtis_rd( const std::vector<double> &dbh,
                          const std::vector<double> &expansion,
                          const bool imperial );                                                                               

double compute_reineke_sdi( const std::vector<double> &dbh,
                            const std::vector<double> &expansion,
                            const bool imperial );  

double compute_ccf( const std::vector<double> &crown_width,
                    const std::vector<double> &expansion,
                    const bool imperial ); 

double compute_R( const std::vector<double> &x, 
                  const std::vector<double> &y,
                  double plotarea,
                  const std::vector<double> &poly_x,
                  const std::vector<double> &poly_y );

double compute_R( const std::vector<double> &x, 
                  const std::vector<double> &y,
                  double plotarea,                  
                  const Point &plotCenter,
                  const double plotRadius );   

std::vector<double> Ripley_Edge_Correction( const std::vector<Point> &trees, 
                                            const double radius,
                                            const std::vector<Point> &plot );                                 

std::vector<double> compute_Hegyi( const std::vector<double> &x, 
                                   const std::vector<double> &y,
                                   const std::vector<double> &dbh,
                                   const std::vector<Point> &plot,
                                   const bool imperial_units );

std::vector<double> compute_Arney_CSI( const std::vector<double> &x,
                                       const std::vector<double> &y,
                                       const std::vector<double> &dbh,
                                       const std::vector<double> &mcw );                                   

MCWPARMS &get_mcw_parms( int  spcd, 
                         bool allow_donor, 
                         bool allow_none );

std::vector<double> compute_mcw( const std::vector<int> &fia,
                                 const std::vector<double> &dbh,
                                 const bool imperial_units,
                                 const bool allow_donor,
                                 const bool allow_none,
                                 const bool hold_domain,
                                 const bool apply_ceiling );

std::vector<double> compute_lcw( const std::vector<int> &fia,
                                 const std::vector<double> &dbh,
                                 const std::vector<double> &cr,
                                 const bool imperial_units,
                                 const bool allow_donor,
                                 const bool allow_none,
                                 const bool hold_domain,
                                 const bool apply_ceiling );

std::vector<double> height_dbh_fit( const std::vector<double> &dbh,
                                    const std::vector<double> &height,
                                    const double bh );

std::vector<double> height_dbh_predict( const std::vector<double> &parameters,
                                        const std::vector<double> &dbh,
                                        const double bh );

std::vector<double> compute_glover_hool( const std::vector<double> &dbh,
                                         const std::vector<double> &expansion,
                                         const bool use_arithmetic,
                                         const bool imperial ); 

std::vector<double> compute_apa( const std::vector<double> &x,
                                 const std::vector<double> &y,
                                 const std::vector<double> &dbh,
                                 const std::array<Point,2> &plot_corners,
                                 const bool weighted ); 

std::vector<std::vector<Point>> get_voronoi_polygons( 
                                    const std::vector<double> &x,
                                    const std::vector<double> &y,
                                    const std::vector<double> &dbh,
                                    const std::array<Point,2> &plot_corners,
                                    const bool weighted );

#endif