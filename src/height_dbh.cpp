
#include <vector>
#include <cmath>
#include <algorithm> 
#include <numeric>
#include <functional> 

// Function to perform linear regression
std::tuple<double,double,double> _linearRegression( const std::vector<double>& x, 
                                                    const std::vector<double>& y ) 
{
    double intercept = 0.0;
    double slope = 0.0;
    size_t n = x.size();

    if( n != y.size() ) 
        return std::tuple(0.0,0.0,0.0);

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_x2 = 0.0;
    double sum_xy = 0.0;
    for( size_t i = 0; i < n; ++i )
    {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    slope = ( n*sum_xy - sum_x*sum_y ) / ( n*sum_x2 - sum_x*sum_x );
    intercept = ( sum_y - slope*sum_x ) / n;

    double sse = 0.0;
    for( size_t i = 0; i < x.size(); ++i ) 
    {
        double est = y[i] - (intercept + slope*x[i]);
        sse += est * est;
    }

    return std::tuple( intercept, slope, sse );
}


double _hd_function( double dbh, double bh, double a, double b, double c ) 
{
    return bh + std::exp( a + b * std::pow( dbh, c ) ); 
}

double rmse( const std::vector<double> &dbh, 
             const std::vector<double> &height, 
             const double bh, 
             double a, double b, double c )
{
    double sse = 0.0;
    for( size_t i = 0; i < dbh.size(); ++i )
    {
        double est = _hd_function( dbh[i], bh, a, b, c );
        sse += (height[i] - est) * (height[i] - est);
    }

    return std::sqrt(sse/dbh.size());
}

std::vector<double> height_dbh_fit( const std::vector<double> &dbh,
                                    const std::vector<double> &height,
                                    const double bh )
{
    if( dbh.size() != height.size() )
        return std::vector<double>{};

    std::vector<double> htadbh( height.size(), 0.0 );
    std::vector<double> dbhc( dbh.size(), 0.0 );
    std::vector<double> parameters(3,0.0);
    double a, b;
    double sse = 0.0;
    double low_rmse = __DBL_MAX__;
   
    std::transform( height.begin(), height.end(), htadbh.begin(), [](double h){ return std::log(h-4.5); });
  
    for( double c = -1.0; c < -0.1; c += 0.01 )
    {
        std::transform( dbh.begin(), dbh.end(), dbhc.begin(), [c](double z){ return std::pow( z, c); });
        std::tie( a, b, sse ) = _linearRegression( dbhc, htadbh );
        double cur_rmse = rmse( dbh, height, bh, a, b, c );

        if( cur_rmse < low_rmse ) 
        {
            low_rmse = cur_rmse;
            parameters = { a, b, c };
        }
    }

    return parameters;
}


std::vector<double> height_dbh_predict( const std::vector<double> &parameters,
                                        const std::vector<double> &dbh,
                                        const double bh )
{
    if( parameters.size() == 0 || dbh.size() == 0 || parameters.size() != 3 )
        return std::vector<double>{};

    std::vector<double> height_hat;

    for( auto &d : dbh )
        height_hat.emplace_back( _hd_function( d, bh, parameters[0], parameters[1], parameters[2] ));

    return height_hat;
}
