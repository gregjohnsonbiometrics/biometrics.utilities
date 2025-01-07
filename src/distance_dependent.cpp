#include "utilities.hpp"
#include <iostream>
#include <cmath>
#include <limits>

struct Point {
    double x, y;
};

// Function to calculate Euclidean distance between two points
double _distance(const Point& p1, const Point& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

std::vector<double> _findNearestNeighborDistance( const std::vector<double> &x, 
                                                  const std::vector<double> &y ) 
{
    int n = x.size();
    std::vector<double> nn_distance( n, 0.0 ); 
    std::vector<size_t>indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    // Sort points by x-coordinate to improve search efficiency
    indices = sort_indices( x );

    for( int i = 0; i < n; ++i ) 
    {
        int k = indices[i];
        double minDist = std::numeric_limits<double>::max();

        // Only compare with nearby points in sorted order to reduce checks
        for( int j = i + 1; j < n && (x[indices[j]] - x[k]) < minDist; ++j ) 
        {
            double dist = _distance( Point(x[k],y[k]), Point(x[indices[j]],y[indices[j]]) );
            if (dist < minDist) 
                minDist = dist;
        }

        // Check points on the left side in sorted order
        for( int j = i - 1; j >= 0 && (x[k] - x[indices[j]]) < minDist; --j ) 
        {
            double dist = _distance( Point(x[k],y[k]), Point(x[indices[j]],y[indices[j]]) );
            if (dist < minDist) 
                minDist = dist;
        }

        nn_distance[k] = minDist;
    }

    return nn_distance;
}


// Function to compute the distance to the nearest edge
double _distanceToEdge( const double x, 
                        const double y, 
                        const double width, 
                        const double height ) 
{
    double dx = std::min(x, width - x);
    double dy = std::min(y, height - y);
    return std::min(dx, dy);
}

// Function to compute Donnelly's edge correction
std::vector<double> _DonnellyCorrection(const std::vector<double> &x, 
                                        const std::vector<double> &y, 
                                        double width, double height, double radius ) 
{
    auto n = x.size();
    std::vector<double> corrections(n, 1.0);

    for( size_t i = 0; i < n; ++i ) 
    {
        double d = _distanceToEdge( x[i], y[i], width, height );
        if (d < radius) {
            corrections[i] = (PI * radius * radius) / (PI * radius * radius - 2 * radius * (radius - d));
        }
    }

    return corrections;
}


double compute_R( const std::vector<double> &x, 
                  const std::vector<double> &y,
                  const double plot_area,
                  const double ulx,
                  const double uly,
                  const double x_width,
                  const double y_width )
{
    double n = x.size();
    auto distances = _findNearestNeighborDistance( x, y );
    auto average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n;

    // if ulx and uly are provided ( ulx and uly != -1 ), then compute area, and point offsets for Donnelly correction
    if( ulx != -1.0 || uly != -1.0 || x_width > 0.0 || y_width > 0.0 )
    {
        // compute x,y offets from upper left corner
        std::vector<double> x_o( n );
        std::vector<double> y_o( n );
        std::transform( x.begin(), x.end(), x_o.begin(), [ulx](double p){return p - ulx;} );
        std::transform( y.begin(), y.end(), y_o.begin(), [uly](double p){return p - uly;} );

        // compute Donnelly bias correction
        auto corrections = _DonnellyCorrection( x_o, y_o, x_width, y_width, average_distance );

        // adjust distances for bias correction
        std::transform( distances.begin(), distances.end(), corrections.begin(), 
                        distances.begin(), [](double a, double b){ return a*b; } );

        // recompute average distance with bias correction
        average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n; 
    }

     auto R = average_distance / (std::sqrt(plot_area/n)/2.0);

    return R;
}

std::vector<double> compute_Hegyi( const std::vector<double> &x, 
                                   const std::vector<double> &y,
                                   const std::vector<double> &dbh,
                                   const bool imperial_units )
{
    size_t n = x.size();

    std::vector<double> h( n, 0.0 );

    for( size_t i = 0; i < n; ++i )
    {
        if( dbh[i] > 0.0 )
        {
            auto di2 = dbh[i]*dbh[i];

            auto xyi = Point(x[i], y[i]);

            for( size_t j = 0; j < n; ++j )
            {
                if( j != i )
                {
                    auto d = _distance( xyi, Point( x[j], y[j] ) );
                    if( imperial_units ) d *= 0.3048;
                    if( d <= 6.0 && d > 0.0 )
                        h[i] += ( (dbh[j]*dbh[j]) / di2 ) / d;
                }
            }
        } else {
            h[i] = NAN;
        }
    }

    return h;
}



std::vector<double> compute_Arney_CSI( const std::vector<double> &x,
                                       const std::vector<double> &y,
                                       const std::vector<double> &dbh,
                                       const std::vector<double> &mcw )
{
    auto n = dbh.size();
    std::vector<double> csi( n, 0.0 );
    std::vector<double> crown_area( n, 0.0 );

    // compute maximum crown area for each tree
    std::transform( mcw.begin(), mcw.end(), crown_area.begin(), [](double mcw_t){ return PI*mcw_t*mcw_t/4.0; } );
   
    for( size_t i = 0; i < n; ++i )
    {
        // compute distances to current tree i
        std::vector<double> distances( n, 0.0 );
        for( size_t j = 0; j < n; ++j )
            distances[j] = _distance( Point(x[i],y[i]), Point(x[j],y[j]) );

        std::vector<double> area_overlap(n,0.0);

        // compute crown overlap
        for( size_t j = 0; j < n; ++j )
        {
            double r1 = 0.0;
            double r2 = 0.0;

            // arrange the largest and smallest tree for AO and CSI
            if( dbh[i] >= dbh[j] )
            {
                // largest tree is subject tree
                r1 = mcw[i]/2.0;

                // smaller tree is competitor
                r2 = mcw[j]/2.0;
            } else {
                // the largest tree is the competitor
                r1 = mcw[j]/2.0;

                // smaller tree is the subject tree
                r2 = mcw[i]/2.0;
            }

            // if this is the subject tree
            if( i == j )
            {
                // overlap area is 0.0
                csi[i] += 100.0;
                continue;
            }

            // check if trees actually overlap
            if( r1 + r2 < distances[j] )
            {
                // no overlap, overlap area is 0.0
                continue;
            }

            // r1 is largest and totally overlaps r2
            if( r1 > r2 && r1 > r2 + distances[j] )
            {
                // arrangement 4: smaller tree is completely overlapped by the larger tree
                // (Equation 4 Appendix 1 page 5)
                area_overlap[j] = PI * r2 * r2;
                csi[i] += 100.0*area_overlap[j]/(PI*r1*r1);
                continue;
            }

            auto s = ( r1 + r2 + distances[j] ) / 2.0;
            auto c = ( 2.0 / distances[j] ) * std::sqrt( s*(s-r1)*(s-r2)*(s-distances[j]) );
            auto x1 = std::sqrt( r1*r1 - c*c );

            /////////// not sure why we have to round here
            //DIST = round(DIST,0)
            //x1 = round(x1,0)
            
            // arrangement 1: distance > x1 (Equation 1 Appendix 1 page 5)
            if( distances[j] > x1 )
            {
                area_overlap[j] = r1*r1*std::asin(c/r1) + r2*r2 * std::asin(c/r2) - distances[j]*c;
                csi[i] += 100.0*area_overlap[j]/(PI*r1*r1);
                continue;
            }

            // arrangement 2: DIST = x1 (Equation 2 Appendix 1 page 5)
            if( distances[j] == x1 )
            {
                area_overlap[j] = (0.5*PI*r2*r2) + (r1*r1*std::asin(r2/r1)) - (distances[j]*r2);
                csi[i] += 100.0*area_overlap[j]/(PI*r1*r1);
                continue;
            }

            if( distances[j] < x1 )
            {
                // arrangement 3: DIST < x1 but competitor is not completely overlapped
                // (Equation 3 Appendix 1 page 5)
                auto x2 = x1 - distances[j];
                area_overlap[j] = PI*r2*r2 - r2*r2*std::asin(c/r2) + x2*c + r1*r1 * std::asin(c/r1) - x1*c;
                csi[i] += 100.0*area_overlap[j]/(PI*r1*r1);
                continue;
            }
        }

        if( i == 0 )
        {
            for( size_t j = 0; j < n; ++j )
            std::cout << area_overlap[j] << "\t" << csi[j] << "\n";
        }
    }

    return csi;
}
