
#include <random>       // std::default_random_engine, std::uniform_real_distribution

#include <math.h>
#include <exception>
#include <stdexcept>

#include <iostream>

#include "utilities.hpp"


double compute_ba( const std::vector<double> &dbh, const std::vector<double> &expansion_factor, const bool imperial )
{
    if( dbh.size() == 0 || dbh.size() != expansion_factor.size() )
        return NAN;

    double _ba = 0.0;
    auto k = imperial ? 0.005454154 : 0.000025*PI;

    // traverse the vector computing and storing bal
    for( auto i : sort_indices(dbh ) )
        _ba += dbh[i] * dbh[i] * expansion_factor[i] * k;

    return _ba;
}


std::vector<double> compute_bal( const std::vector<double> &dbh, const std::vector<double> &expansion_factor, const bool imperial )
{
    if( dbh.size() == 0 || dbh.size() != expansion_factor.size() )
        return std::vector<double>{};

    double _bal = 0.0;
    double pending = 0.0;
    std::vector<double>bal(dbh.size(),0.0);
    double last_dbh = 9999.0;
    auto k = imperial ? 0.005454154 : 0.000025*PI;

    // traverse the vector computing and storing bal
    for( auto i : sort_indices(dbh ) )
    {
        double tree_ba = dbh[i] * dbh[i] * expansion_factor[i] * k;
 
        if( dbh[i] < last_dbh ) 
        {
            bal[i] = _bal;
            pending = _bal;
            _bal += tree_ba;
            last_dbh = dbh[i];
        } else if( dbh[i] == last_dbh ) {
            bal[i] = pending;
            _bal += tree_ba;
        } else {
            throw std::invalid_argument("dbh vector sort error in compute_bal\n");
        }
    }

    return bal;
}


std::vector<double> compute_ccfl( const std::vector<double> &dbh, const std::vector<double> &mcw, const std::vector<double> &expansion, const bool imperial=true )
{
    if( dbh.size() == 0 || dbh.size() != mcw.size() || dbh.size() != expansion.size() )
        return std::vector<double>{};

    double _ccfl = 0.0;
    double pending = 0.0;
    std::vector<double>ccfl(dbh.size(),0.0);
    double last_dbh = 9999.0;

    double area = imperial ? 43560.0 : 10000.0;

    // traverse the vector computing and storing bal
    for( auto i : sort_indices(dbh ) )
    {
        auto mca = 100.0*((PI*(mcw[i]*mcw[i]/4.0))/area)*expansion[i];
        if( dbh[i] < last_dbh ) 
        {
            ccfl[i] = _ccfl;
            pending = _ccfl;
            _ccfl += mca;
            last_dbh = dbh[i];
        } else if( dbh[i] == last_dbh ) {
            ccfl[i] = pending;
            _ccfl += mca;
        } else {
            throw std::invalid_argument("dbh vector sort error in compute_ccfl\n");
        }
    }

    return ccfl;
}

// compute dominant height
// method 0: average height of the dominant_cohort_size trees by decreasing dbh
//        1: average height of the dominant_cohort_size trees by decreasing height
//        2: Lorey height (height of the tree of average basal area) 
double compute_dominant_height( const std::vector<double> &height,
                                const std::vector<double> &dbh,
                                const std::vector<double> &expansion,
                                const int dominant_cohort_size,
                                const int method )
{
    if( height.size() == 0 || dbh.size() != expansion.size() || dbh.size() != height.size() )
        return NAN;

    double dominant_height = 0.0;

    if( method == 2 )
    {
        // compute lorey height
        double sum_ht = 0.0;
        double sum_ba = 0.0;
        for( size_t i = 0; i < height.size(); i++ )
        {
            sum_ht += height[i] * dbh[i] * dbh[i] * expansion[i];
            sum_ba += dbh[i] * dbh[i] * expansion[i];
        }

        dominant_height = (sum_ba > 0.0) ? sum_ht / sum_ba : 0.0;
    } else if( method < 2 ) {

        if( dominant_cohort_size <= 0.0 )
            return NAN;
            
        double accumulated_expansion = 0.0;
        double sum_height = 0.0;

        for( auto i : sort_indices( method == 0 ? dbh : height ) )
        {  
            if( accumulated_expansion + expansion[i] <= dominant_cohort_size )
            {
                sum_height += height[i] * expansion[i];
                accumulated_expansion += expansion[i];
            } else if( accumulated_expansion + expansion[i] > dominant_cohort_size ) {
                double n_trees = (double)dominant_cohort_size - accumulated_expansion;
                sum_height += height[i] * n_trees;
                accumulated_expansion = dominant_cohort_size;
                break;
            }
        }

        dominant_height = (accumulated_expansion > 0.0) ? sum_height / accumulated_expansion : 0.0;        
    } else {
        throw std::invalid_argument("invalid method in compute_dominant_height\n");
    }

    return dominant_height;
}

// compute dbh of dominant height cohort
// method 0: average dbh of the dominant_cohort_size trees by decreasing dbh
//        1: average dbh of the dominant_cohort_size trees by decreasing height
//        2: Lorey dbh (quadratic mean diameter)
double compute_dominant_dbh( const std::vector<double> &height,
                             const std::vector<double> &dbh,
                             const std::vector<double> &expansion,
                             const int dominant_cohort_size,
                             const int method )
{
    if( height.size() == 0 || dbh.size() != expansion.size() || dbh.size() != height.size() )
        return NAN;

    double dominant_dbh = 0.0;

    if( method == 2 )
    {
        // compute lorey height
        double sum_expansion = 0.0;
        double sum_ba = 0.0;
        for( size_t i = 0; i < height.size(); i++ )
        {
            sum_expansion += expansion[i];
            sum_ba += dbh[i] * dbh[i] * expansion[i];
        }

        dominant_dbh = (sum_expansion > 0.0) ? std::sqrt( sum_ba / sum_expansion ) : 0.0;

    } else if( method < 2 ) {

        if( dominant_cohort_size <= 0.0 )
            return NAN;
            
        double accumulated_expansion = 0.0;
        double sum_dbh = 0.0;

        for( auto i : sort_indices( method == 0 ? dbh : height ) )
        {  
            if( accumulated_expansion + expansion[i] <= dominant_cohort_size )
            {
                sum_dbh += dbh[i] * expansion[i];
                accumulated_expansion += expansion[i];
            } else if( accumulated_expansion + expansion[i] > dominant_cohort_size ) {
                double n_trees = (double)dominant_cohort_size - accumulated_expansion;
                sum_dbh += dbh[i] * n_trees;
                accumulated_expansion = dominant_cohort_size;
                break;
            }
        }

        dominant_dbh = (accumulated_expansion > 0.0) ? sum_dbh / accumulated_expansion : 0.0;        
    } else {
        throw std::invalid_argument("invalid method in compute_dominant_dbh\n");
    }

    return dominant_dbh;
}


double compute_qmd( const std::vector<double> &dbh,
                    const std::vector<double> &expansion )
{
    if( dbh.size() == 0 || dbh.size() != expansion.size() )
        return NAN;

    // compute quadratic mean diameter
    double stocking = 0.0;
    double qmd = 0.0;
    for( size_t i = 0; i < dbh.size(); i++ )
    {
        stocking += expansion[i];
        qmd += dbh[i] * dbh[i] * expansion[i];
    }

    if( stocking > 0.0 )
    {
        qmd = std::sqrt( qmd/stocking );
    } else {
        throw std::invalid_argument("stocking is 0.\n");    
    }

    return qmd;
} 

double compute_relative_spacing( const std::vector<double> &expansion,
                                 const double dominant_height,
                                 const bool imperial )
{
    if( expansion.size() == 0 || dominant_height <= 0.0 )
        return NAN;

    double area = (imperial ? 43560.0 : 10000.0 );
    auto stocking = std::accumulate( expansion.begin(), expansion.end(), 0.0 );

    double rs = 0.0;
    
    if( stocking > 0.0 && dominant_height > 0.0 )
        rs = std::sqrt( area / stocking ) / dominant_height;

    return rs;
}

double compute_curtis_rd( const std::vector<double> &dbh,
                          const std::vector<double> &expansion,
                          const bool imperial )
{
    if( expansion.size() == 0 || dbh.size() != expansion.size() )
        return NAN;

    double k = imperial ? 0.005454154 : 0.00007853975;

    // compute stocking and quadratic mean diameter
    double stocking = 0.0;
    double qmd = 0.0;
    for( size_t i = 0; i < dbh.size(); i++ )
    {
        stocking += expansion[i];
        qmd += dbh[i] * dbh[i] * expansion[i];
    }

    double curtis_rd = 0.0;
    if( stocking > 0.0 )
    {
        double ba = qmd * k;
        qmd = std::sqrt( qmd/stocking );
        curtis_rd = ba / std::sqrt( qmd );
    } 

    return curtis_rd;
} 

double compute_reineke_sdi( const std::vector<double> &dbh,
                            const std::vector<double> &expansion,
                            const bool imperial )
{
    if( expansion.size() == 0 || dbh.size() != expansion.size() )
        return NAN;

    double k = imperial ? 10.0 : 25.4;

    // compute stocking and quadratic mean diameter
    double stocking = 0.0;
    double qmd = 0.0;
    for( size_t i = 0; i < dbh.size(); i++ )
    {
        stocking += expansion[i];
        qmd += dbh[i] * dbh[i] * expansion[i];
    }

    double sdi = 0.0;
    if( stocking > 0.0 )
    {
        qmd = std::sqrt( qmd/stocking );
        sdi = stocking * std::pow(qmd/k,1.605);
    } 

    return sdi;
} 

double compute_ccf( const std::vector<double> &crown_width,
                    const std::vector<double> &expansion,
                    const bool imperial )
{
    if( expansion.size() == 0 || crown_width.size() != expansion.size() )
        return NAN;

    const double area = imperial ? 43560.0 : 10000.0;

    double ca = 0.0;
    for( size_t i = 0; i < crown_width.size(); i++ )
        ca += (crown_width[i] * crown_width[i])/4.0 * PI * expansion[i];

    double ccf = ca / area * 100.0;

    return ccf;
}


std::vector<double> compute_glover_hool( const std::vector<double> &dbh,
                                         const std::vector<double> &expansion,
                                         const bool use_arithmetic, 
                                         const bool imperial )
{
    auto k = imperial ? 0.005454154 : 0.000025*PI;
    std::vector<double> G;

    if( expansion.size() == 0 || dbh.size() != expansion.size() )
        return std::vector<double>{};

    double mean = 0.0;
    if( use_arithmetic )
    {
        double n = 0.0;
        for( size_t i = 0; i < dbh.size(); ++i )
        {
            mean += dbh[i] * expansion[i];
            n += expansion[i];
        }

        if( n > 0.0 )
            mean /= n;
        else
            return std::vector<double>{};

    } else {
        mean = compute_qmd( dbh, expansion );
    }
    
    if( mean == 0.0 )
        return std::vector<double>{};

    for( auto &d : dbh )
        G.push_back( d*d*k / mean );

    return G;
}


// upper-crown exponents (x, y). ORGANON rows are the refit; SW and HW are the
// shipped forestCI defaults (gothic and ellipsoid).

const std::unordered_map<int, std::pair<double, double>> known_shape = {{
    { 202, { 1.099, 1.258 } },
    { 263, { 0.710, 1.678 } },
    { 352, { 0.854, 1.957 } }
    }};

const std::unordered_map<int, std::pair<double, double>> default_shape = {{
    {   0, { 3.000, 1.000 } }, // softwood
    {   1, { 3.000, 3.000 } }  // hardwood
    }};

std::pair<double, double> get_shape( int fia ) 
{
    std::pair<double, double> s;
    
    try {
        s =  known_shape.at( fia );
    } catch( std::exception &e ) {
        s = default_shape.at( fia < 300 ? 0 : 1 );
    }
    
    return s;
}

// widest-point rule: proportion of crown length above crown base at which
// the crown is widest (0 = at crown base). 

double widest_point( int fia,
                     double dbh,
                     double height,
                     double crown_length,
                     bool imperial_units )
{
    double lcl_rel = 0.0;
    
    double to_ft = imperial_units ? 1.0 : 3.280840;
    double to_in = imperial_units ? 1.0 : 0.3937008;
    
    switch( fia ) {
    case 202: //  Douglas-fir: Hann and Mainwaring 2026, Eq 20
        lcl_rel = 1.0 - std::exp( -0.699472314 - 0.548622967 * std::sqrt( dbh * to_in ) +
            0.00962202678 * height * to_ft );
        
        lcl_rel = std::min( std::max( 1.0 - lcl_rel, 0.0 ), 0.999 );
        break;
        
    case 263: // Western Hemlock: Marshall et al. 2003, Eq 6
        lcl_rel = 0.355270;
        break;
        
    case 352: // red alder Hann et al. 2011, Table 6.2
        lcl_rel = 0.63619616 * std::exp( -1.2180562 * std::pow( (1.0 - height * to_ft / 140.0 ), 3.0 ) );
        lcl_rel = std::min( std::max( lcl_rel, 0.0 ), 0.999 );
        break;
        
    default: 
        if( fia < 300 )
            lcl_rel = 0.20; // forestCI softwood default
    else
        lcl_rel = 0.40; // forestCI hardwood default
    }
    
    return lcl_rel;
}


double compute_crown_profile( 
                      const int fia, 
                      const double dbh, 
                      const double height,
                      const double hcb, 
                      const double lcw,
                      const double h, 
                      const bool imperial_units ) 
{
    // compute crown length
    double cl = height - hcb;
    
    // get widest point (relative) in crown ( 0 is base, 1 is tip )
    double w = widest_point( fia, dbh, height, cl, imperial_units );
    
    // height of largest crown width
    double hlcw = hcb + w * cl;
    
    double z = std::min( std::max( (h - hlcw) / (height - hlcw), 0.0 ), 1.0 );
    auto [x, y] = get_shape( fia );
    
    // compute crown width at h
    double cwa = (h < hlcw) ? 0.0 : lcw * std::pow( std::max( 1.0 - std::pow( z, x ), 0.0 ), ( 1.0 / y ) );
    
    return cwa;
}
