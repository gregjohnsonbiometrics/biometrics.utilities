#include "utilities.hpp"

#include <iostream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <ranges>


// Function to find translation values to the origin
std::pair<double,double> _translateToOrigin( const std::vector<Point> &plot ) 
{
    auto min_x = std::min_element( plot.begin(), plot.end(), [](const auto a, const auto b){ return a.x < b.x;} )->x;
    auto min_y = std::min_element( plot.begin(), plot.end(), [](const auto a, const auto b){ return a.y < b.y;} )->y;

    return std::pair( min_x, min_y );
}

// Function to rotate the polygon around the origin
std::vector<Point> _rotatePolygon( const std::vector<Point>& polygon, double angle ) 
{
    std::vector<Point> new_polygon = polygon;

    double radians = angle * PI / 180.0; // Convert angle to radians
    double cosTheta = std::cos(radians);
    double sinTheta = std::sin(radians);
    
    for( auto &point : new_polygon ) 
    {
        double xNew = point.x * cosTheta - point.y * sinTheta;
        double yNew = point.x * sinTheta + point.y * cosTheta;
        point.x = xNew;
        point.y = yNew;
    }

    return new_polygon;
}

struct MinMax {
    double minX;
    double maxX;
    double minY;
    double maxY;
};

MinMax _calculateBoundingBox(const std::vector<Point> &plot,
                             const std::pair<double,double> &trans,
                             double angle ) 
{
    MinMax mm = { 0.0, 0.0, 0.0, 0.0 };
    auto &[x_trans, y_trans] = trans;
    mm.minX = mm.maxX = plot[0].x - x_trans;
    mm.minY = mm.maxY = plot[0].y - y_trans;
    
    double radians = angle * PI / 180.0;
    double cosTheta = std::cos(radians);
    double sinTheta = std::sin(radians);

    for( const auto& point : plot ) 
    {
        double xRot = (point.x - x_trans) * cosTheta - (point.y - y_trans) * sinTheta;
        double yRot = (point.x - x_trans) * sinTheta + (point.y - y_trans) * cosTheta;

        mm.minX = std::min(mm.minX, xRot);
        mm.maxX = std::max(mm.maxX, xRot);
        mm.minY = std::min(mm.minY, yRot);
        mm.maxY = std::max(mm.maxY, yRot);
    }

    return mm;
}


// Function to find the optimal rotation angle to align the bounding box with axes
double _findOptimalRotationAngle( const std::vector<Point> &plot,
                                  const std::pair<double,double> &trans ) 
{
    double bestAngle = 0.0;
    double minArea = 1.0e9;
    
    for( double angle = 0; angle < 180.0; angle += 1.0 ) 
    {
        MinMax mm = _calculateBoundingBox( plot, trans, angle );
        double area = ( mm.maxX - mm.minX ) * ( mm.maxY - mm.minY) ;
        
        if( area < minArea ) 
        {
            minArea = area;
            bestAngle = angle;
        }
    }
    
    return bestAngle;
}

double _calculatePartArea( double r, double w, double m, double n1, double n2 )
{
    double r2 = r * r;

    if( m < 0.01 )
        return 0.0;

    double theta = std::atan(n1 / m) + std::atan(n2 / m);

    if (m >= r )
    {
        // case a
        return r * r * theta / 2;
    }

    double m2 = m * m;

    if( n1 * n1 + m2 >= r2 )
    {
        if( n2 * n2 + m2 >= r2 )
        {
            // case b
            return r2 * (theta - 2 * std::acos(m / r)) / 2.0  + m * std::sqrt(r2 - m2);
        }

        // case d
        double d = std::sqrt(r2 - m2) + n2;
        return d * m / 2.0 + r2 * (theta - std::atan(n2 / m) - std::acos(m / r)) / 2.0;
    }

    if( n2 * n2 + m2 >= r2 )
    {
        // case c
        double d = std::sqrt(r2 - m2) + n1;
        return d * m / 2.0 + r2 * (theta - std::atan(n1 / m) - std::acos(m / r)) / 2.0;
    }

    // case e
    return m * w / 2;
}


double _calculateArea( Point p, 
                       double r, 
                       double w, 
                       double h )
{
    if (p.x < 0 || p.y < 0 || p.x > w || p.y > h )
    {
        throw("The circle's center must be contained by the rectangle.");
    }

    return _calculatePartArea(r, w, p.y, p.x, w - p.x)
         + _calculatePartArea(r, w, p.x, p.y, h - p.y)
         + _calculatePartArea(r, w, h - p.y, p.x, w - p.x)
         + _calculatePartArea(r, w, w - p.x, p.y, h - p.y);
}

// Computes Ripley's (1977) isotropic edge correction weight
std::vector<double> Ripley_Edge_Correction( const std::vector<Point> &trees,    // tree locations
                                            const double radius,                // influence radius to test for edge overlap
                                            const std::vector<Point> &plot )    // plot corner coordinates
{
    std::vector<double> weights( trees.size(), 1.0 );

    double circle_area = radius * radius * PI; 
    auto trans = _translateToOrigin( plot );

    double optimalAngle = _findOptimalRotationAngle( plot, trans );

    auto rotated_plot = _rotatePolygon( plot, optimalAngle );

    auto width =  std::max_element( rotated_plot.begin(), rotated_plot.end(), [](const auto a, const auto b){ return a.x < b.x;} )->x;
    auto height = std::max_element( rotated_plot.begin(), rotated_plot.end(), [](const auto a, const auto b){ return a.y < b.y;} )->y;

    for( size_t i = 0; i < trees.size(); ++i )
        weights[i] =  circle_area / _calculateArea( trees[i], radius, width, height ); 

    return weights;
}
