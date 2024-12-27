#include "utilities.hpp"
#include <iostream>
#include <cmath>
#include <limits>

struct Point {
    double x, y;
};

// Function to calculate Euclidean distance between two points
double distance(const Point& p1, const Point& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

std::vector<double> findNearestNeighborDistance( std::vector<double> x, 
                                                 std::vector<double> y ) 
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
            double dist = distance( Point(x[k],y[k]), Point(x[indices[j]],y[indices[j]]) );
            if (dist < minDist) 
                minDist = dist;
        }

        // Check points on the left side in sorted order
        for( int j = i - 1; j >= 0 && (x[k] - x[indices[j]]) < minDist; --j ) 
        {
            double dist = distance( Point(x[k],y[k]), Point(x[indices[j]],y[indices[j]]) );
            if (dist < minDist) 
                minDist = dist;
        }

        nn_distance[k] = minDist;
    }

    return nn_distance;
}


double compute_R( const std::vector<double> x, 
                  const std::vector<double> y,
                  const double plot_area )
{
    auto distances = findNearestNeighborDistance( x, y );

    double n = x.size();

    auto average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n;

    auto R = average_distance / (std::sqrt(plot_area/n)/2.0);

    return R;
}
