#include "utilities.hpp"
#include "FortuneAlgorithm.h"
#include <iostream>
#include <limits>

constexpr double OFFSET = 1.0f;

// Function to calculate Euclidean distance between two points
double _distance(const Point& p1, const Point& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

// compute area of a polygon using Shoelace algorithm
double _computePolygonArea( const std::vector<Point>& vertices ) 
{
    int n = vertices.size();
    if (n < 3) 
       throw( "A polygon must have at least 3 vertices.\n" );

    double area = 0.0;

    // Apply Shoelace formula
    for (int i = 0; i < n; ++i) 
    {
        int next = (i + 1) % n; // Ensure it wraps around to the first vertex
        area += (vertices[i].x * vertices[next].y) - (vertices[next].x * vertices[i].y);
    }

    return std::abs(area) / 2.0;
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
            double dist = _distance( {x[k],y[k]}, {x[indices[j]],y[indices[j]]} );
            if (dist < minDist) 
                minDist = dist;
        }

        // Check points on the left side in sorted order
        for( int j = i - 1; j >= 0 && (x[k] - x[indices[j]]) < minDist; --j ) 
        {
            double dist = _distance( {x[k],y[k]}, {x[indices[j]],y[indices[j]]} );
            if (dist < minDist) 
                minDist = dist;
        }

        nn_distance[k] = minDist;
    }

    return nn_distance;
}

// Function to compute the minimum distance from a point to the edges of the plot
double _minDistanceToPolygonEdge( const Point& p, 
                                  const std::vector<Point>& polygon ) 
{
    size_t n = polygon.size();
    double minDist = std::numeric_limits<double>::max();

    for( size_t i = 0; i < n; ++i ) 
    {
        Point p1 = polygon[i];
        Point p2 = polygon[(i + 1) % n];

        // Compute the distance from point p to the edge (p1, p2)
        double A = p.y - p1.y;
        double B = p.x - p1.x;
        double C = p2.x - p1.x;
        double D = p2.y - p1.y;

        double dot = A * C + B * D;
        double len_sq = C * C + D * D;
        double param = (len_sq != 0) ? (dot / len_sq) : -1;

        Point nearest;
        if (param < 0) {
            nearest = p1;
        } else if (param > 1) {
            nearest = p2;
        } else {
            nearest = {p1.x + param * C, p1.y + param * D};
        }

        double dist = _distance(p, nearest);
        minDist = std::min(minDist, dist);
    }

    return minDist;
}


// Function to compute the distance from a point to the edges of a circular plot
double _DistanceToCircleEdge( const Point &p, 
                              const Point &center,
                              const double radius ) 
{
    double dist = radius - _distance( p, center );

    return dist;
}


// Function to compute Donnelly's edge correction -- polygonal plot version
std::vector<double> _DonnellyCorrection(const std::vector<double> &x, 
                                        const std::vector<double> &y, 
                                        const std::vector<Point> &plotPolygon,
                                        double radius ) 
{
    auto n = x.size();
    std::vector<double> corrections(n, 1.0);

    for( size_t i = 0; i < n; ++i ) 
    {
        double d = _minDistanceToPolygonEdge( {x[i],y[i]}, plotPolygon );
        
        if (d < radius) {
            corrections[i] = (PI * radius * radius) / (PI * radius * radius - 2 * radius * (radius - d));
        }
    }

    return corrections;
}


// Function to compute Donnelly's edge correction -- circle plot version
std::vector<double> _DonnellyCorrection(const std::vector<double> &x, 
                                        const std::vector<double> &y, 
                                        const Point &plotCenter,
                                        const double plotRadius,
                                        double radius ) 
{
    auto n = x.size();
    std::vector<double> corrections(n, 1.0);

    for( size_t i = 0; i < n; ++i ) 
    {
        double d = _DistanceToCircleEdge( {x[i],y[i]}, plotCenter, plotRadius);
        
        if (d < radius) {
            corrections[i] = (PI * radius * radius) / (PI * radius * radius - 2 * radius * (radius - d));
        }
    }

    return corrections;
}

// polygonal plot version
double compute_R( const std::vector<double> &x, 
                  const std::vector<double> &y,
                  double plotarea,
                  const std::vector<double> &poly_x,
                  const std::vector<double> &poly_y )
{
    double n = x.size();

    // if no area and no polygon information, cannot compute R
    if( plotarea <= 0.0 && (poly_x.size() == 0 || poly_y.size() == 0 ) )
        return NAN;

    // if mismatch on vector sizes, return NAN
    if( n == 0 || y.size() == 0 || n != y.size() )
        return NAN;

    auto distances = _findNearestNeighborDistance( x, y );
    auto average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n;
    std::vector<Point> plotPolygon( poly_x.size() );

    if( poly_x.size() > 0 && poly_y.size() > 0 )
    {
        // build polygon vector
        for( size_t i = 0; i < poly_x.size(); ++i )
            plotPolygon[i] = { poly_x[i], poly_y[i] };

        // recompute plot area with supplied polygon
        plotarea = _computePolygonArea( plotPolygon );

        // compute Donnelly bias correction
        auto corrections = _DonnellyCorrection( x, y, plotPolygon, average_distance );

        // adjust distances for bias correction
        std::transform( distances.begin(), distances.end(), corrections.begin(), 
                        distances.begin(), [](double a, double b){ return a*b; } );
    } 

    // recompute average distance with bias correction
    average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n; 

    double R = average_distance / (std::sqrt(plotarea/n)/2.0);

    return R;
}

// circle plot version
double compute_R( const std::vector<double> &x, 
                  const std::vector<double> &y,
                  double plotarea,
                  const Point &plotCenter,
                  const double plotRadius )
{
    double n = x.size();

    // if mismatch on vector sizes, return NAN
    if( n == 0 || y.size() == 0 || n != y.size() )
        return NAN;

    auto distances = _findNearestNeighborDistance( x, y );
    auto average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n;

    if( plotRadius > 0 )
    {
        // recompute plot area with circle dimensions
        plotarea = PI * plotRadius*plotRadius;

        // compute Donnelly bias correction
        auto corrections = _DonnellyCorrection( x, y, plotCenter, plotRadius, average_distance );

        // adjust distances for bias correction
        std::transform( distances.begin(), distances.end(), corrections.begin(), 
                        distances.begin(), [](double a, double b){ return a*b; } );
    }

    // compute average distance 
    average_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / n; 

    double R = average_distance / (std::sqrt(plotarea/n)/2.0);

    return R;
}


std::vector<double> compute_Hegyi( const std::vector<double> &x, 
                                   const std::vector<double> &y,
                                   const std::vector<double> &dbh,
                                   const std::vector<Point> &plot,
                                   const bool imperial_units )
{
    size_t n = x.size();

    // if no tree location data or if mismatch on vector sizes, return empty vector
    if( n == 0 || y.size() == 0  || n != y.size() || dbh.size() != n )
        return std::vector<double>{};

    constexpr double radius = 6.0;
    std::vector<double> h( n, 0.0 );
    std::vector<Point> trees( n );
    std::vector<double> weights( n, 1.0 );

    for( size_t i = 0; i < n; ++i )
        trees[i] = { x[i], y[i] };

    // compute Ripley edge correction weights if a plot boundary polygon is present
    if( plot.size() > 0 )
        weights = Ripley_Edge_Correction( trees, radius, plot );

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
                    if( d <= radius && d > 0.0 )
                        h[i] += weights[i] * ( (dbh[j]*dbh[j]) / di2 ) / d;
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
    // if no tree location data or if mismatch on vector sizes, return empty vector
    auto n = x.size();
    if( n == 0 || y.size() == 0 || y.size() != n || dbh.size() != n || mcw.size() != n )
        return std::vector<double>{};

    std::vector<double> crown_area( n, 0.0 );
    std::vector<double> csi( n, 0.0 );  

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
    }

    return csi;
}


std::vector<std::vector<Point>> _get_voronoi( const std::vector<double> &x,
                                              const std::vector<double> &y,
                                              std::vector<double> weights,
                                              const std::array<Point,2> &plot )
{
    size_t n = x.size();
    std::vector<std::vector<Point>> vresult(n);
    std::vector<Vector2> points(n);

    for( size_t i = 0; i < x.size(); ++i )
        points[i] = { x[i], y[i] };

   // Construct diagram
    FortuneAlgorithm algorithm(points, weights );
    algorithm.construct();

    // Bound the diagram
    algorithm.bound(Box{plot[0].x-0.05, plot[0].y-0.05, plot[1].x+0.05, plot[1].y+0.05}); // Make the bounding box slightly bigger than the intersection box
    VoronoiDiagram diagram = algorithm.getDiagram();

   // Intersect the diagram with a box
    bool valid = diagram.intersect(Box{plot[0].x, plot[0].y, plot[1].x, plot[1].y});
    if (!valid)
        throw std::runtime_error("An error occured in the box intersection algorithm");

    for( std::size_t i = 0; i < diagram.getNbSites(); ++i )
    {
        const VoronoiDiagram::Site* site = diagram.getSite(i);
        Vector2 center = site->point;
        VoronoiDiagram::Face* face = site->face;
        VoronoiDiagram::HalfEdge* halfEdge = face->outerComponent;
        if (halfEdge == nullptr)
            continue;
        while (halfEdge->prev != nullptr)
        {
            halfEdge = halfEdge->prev;
            if (halfEdge == face->outerComponent)
                break;
        }
        VoronoiDiagram::HalfEdge* start = halfEdge;
        while (halfEdge != nullptr)
        {
            if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
            {
                Vector2 origin = (halfEdge->origin->point - center) * OFFSET + center;
                vresult[i].emplace_back( Point{origin.x, origin.y} );
            }
            halfEdge = halfEdge->next;
            if (halfEdge == start)
            {
                Vector2 origin = (halfEdge->origin->point - center) * OFFSET + center;
                vresult[i].emplace_back( Point{origin.x, origin.y} );
                break;
            }
               
        }
    }

    return vresult;
}

std::vector<double> compute_apa( const std::vector<double> &x,
                                 const std::vector<double> &y,
                                 const std::vector<double> &dbh,
                                 const std::array<Point,2> &plot_corners,
                                 const bool weighted )
{
    size_t n = x.size();
    if( n == 0 || y.size() == 0 || y.size() != n || dbh.size() != n || plot_corners.size() != 2 )
        return std::vector<double>{};

    std::vector<double> weights( n, 1.0 );
    if( weighted )
        for( size_t i = 0; i < n; ++i )
            weights[i] = dbh[i];

    // get Voronoi polygons for each tree
    std::vector<std::vector<Point>> apa_poly = _get_voronoi( x, y, weights, plot_corners );

    // compute area of each Voronoi polygon
    std::vector<double> voronoi_area( n, 0.0 );
    for( size_t i = 0; i < n; ++i )
        voronoi_area[i] = _computePolygonArea( apa_poly[i] );

    return voronoi_area;
}


std::vector<std::vector<Point>> get_voronoi_polygons( 
                                    const std::vector<double> &x,
                                    const std::vector<double> &y,
                                    const std::vector<double> &dbh,
                                    const std::array<Point,2> &plot_corners,
                                    const bool weighted )
{
    size_t n = x.size();
    if( n == 0 || y.size() == 0 || y.size() != n || dbh.size() != n || plot_corners.size() != 2 )
        return std::vector<std::vector<Point>>{};

    std::vector<double> weights( n, 1.0 );
    if( weighted )
        for( size_t i = 0; i < n; ++i )
            weights[i] = dbh[i];

    // get Voronoi polygons for each tree
    std::vector<std::vector<Point>> apa_poly = _get_voronoi( x, y, weights, plot_corners );

    return apa_poly;
}