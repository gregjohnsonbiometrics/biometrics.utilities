#ifndef BIOUTIL
#define BIOUTIL

#include <vector>
#include <numeric>      // std::iota, std::accumulate
#include <algorithm>    // std::sort, std::stable_sort

constexpr auto PI = 3.14159265358979323846;

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


std::vector<double> compute_bal( const std::vector<double> &dbh, const std::vector<double> &expansion_factor, const bool imperial );

std::vector<double> compute_ccfl( const std::vector<double> &dbh, const std::vector<double> &mcw, const std::vector<double> &expansion, const bool imperial );

double compute_cch( const double ht,
                    const std::vector<double> &dbh,
                    const std::vector<double> &height,
                    const std::vector<double> &crown_length,
                    const std::vector<double> &dacb,        
                    const std::vector<double> &lcw,
                    const std::vector<double> &expansion,
                    const std::vector<double> &parameters,
                    bool imperial );

double compute_dominant_height( const std::vector<double> height,
                                const std::vector<double> dbh,
                                const std::vector<double> expansion,
                                const int dominant_cohort_size,
                                const int method );  

double compute_relative_spacing( const std::vector<double> expansion,
                                 const double dominant_height,
                                 bool imperial );  

double compute_curtis_rd( const std::vector<double> dbh,
                          const std::vector<double> expansion,
                          bool imperial );                                                                               

double compute_reineke_sdi( const std::vector<double> dbh,
                            const std::vector<double> expansion,
                            bool imperial );  

double compute_ccf( const std::vector<double> crown_width,
                    const std::vector<double> expansion,
                    bool imperial );                                                       
#endif