#include "utilities.hpp"
#include<unordered_map>
#include<string>
#include <math.h>
#include <exception>
#include <stdexcept>

// maximum crown width parameters

std::unordered_map<int, MCWPARMS> mcw_parms = {{
    {  11, { "Pacific silver fir", "Abies", 1, true, 4, 1.5, 0, "Smith" }},                              
    {  12, { "balsam fir", "Abies", 2, false, 1.37, 0.572, 0, "Russell & Weiskittel " }},                
    {  15, { "white fir", "Abies", 1, true, 6.188, 1.0069, 0, "Paine & Hann" }},                         
    {  17, { "grand fir", "Abies", 1, true, 6.188, 1.0069, 0, "Paine & Hann" }},                         
    {  19, { "subalpine fir", "Abies", 1, true, 4, 1.5, 0, "Smith" }},                                   
    {  20, { "California red fir", "Abies", 1, true, 3.0884, 0.7871, 0, "Paine & Hann" }},               
    {  21, { "Shasta red fir", "Abies", 1, true, 3.0884, 0.7871, 0, "Paine & Hann" }},                   
    {  22, { "noble fir", "Abies", 1, true, 3.6883, 0.8627, 0, "Paine & Hann" }},                        
    {  41, { "Port-Orford-cedar", "Chamaecyparis", 1, true, 4, 1.65, 0, "Use RA (Smith)" }},             
    {  42, { "Alaska yellow-cedar", "Chamaecyparis", 1, true, 4, 1.65, 0, "Use RA (Smith)" }},           
    {  73, { "western larch", "Larix", 1, true, 3, 1.37, 0, "Smith" }},                                  
    {  81, { "incense-cedar", "Libocedrus", 1, true, 3.2837, 1.2031, -0.0071858, " Paine & Hann" }},     
    {  93, { "Engelmann spruce", "Picea", 1, true, 4, 1.2, 0, "Smith" }},                                
    {  94, { "white spruce", "Picea", 2, false, 1.5, 0.496, 0, "Russell & Weiskittel " }},               
    {  95, { "black spruce", "Picea", 2, false, 0.535, 0.742, 0, "Russell & Weiskittel " }},             
    {  97, { "red spruce", "Picea", 2, false, 1.8, 0.461, 0, "Russell & Weiskittel " }},                 
    {  98, { "Sitka spruce", "Picea", 1, true, 6.5, 1.8, 0, "Smith" }},                                  
    { 108, { "lodgepole pine", "Pinus", 1, true, 3.2095, 1.2633, -0.0082873, " Paine & Hann" }},        
    { 116, { "Jeffrey pine", "Pinus", 1, true, 5.303, 0.9293, 0, "Paine & Hann" }},                     
    { 117, { "sugar pine", "Pinus", 1, true, 4.6601, 1.0702, 0, "Paine & Hann" }},                      
    { 119, { "western white pine", "Pinus", 1, true, 3.2095, 1.2633, -0.0082873, " Paine & Hann" }},    
    { 122, { "ponderosa pine", "Pinus", 1, true, 3.4835, 1.343, -0.0082544, " Paine & Hann" }},         
    { 129, { "eastern white pine", "Pinus", 2, false, 1.24, 0.585, 0, "Russell & Weiskittel " }},       
    { 131, { "loblolly pine", "Pinus", 1, false, 0.738, 0.245, 0.000809, "Smith & Others (FVS)" }},     
    { 202, { "Douglas-fir", "Pseudotsuga", 1, true, 4.7071, 1.8917, -0.01132, " Modified Arney" }},     
    { 231, { "Pacific yew", "Taxus", 1, true, 4, 1.65, 0, "Use RA (Smith)" }},                          
    { 241, { "northern white-cedar", "Thuja", 2, false, 1.63, 0.436, 0, "Russell & Weiskittel " }},     
    { 242, { "western redcedar", "Thuja", 1, true, 4, 1.65, 0, "Smith" }},                              
    { 261, { "eastern hemlock", "Tsuga", 2, false, 2.44, 0.408, 0, "Russell & Weiskittel " }},          
    { 263, { "western hemlock", "Tsuga", 1, true, 4.5652, 1.4147, 0, "Paine & Hann" }},                 
    { 264, { "mountain hemlock", "Tsuga", 1, true, 6.71, 0.421, 0, "Warbington & Levitan" }},           
    { 312, { "bigleaf maple", "Acer", 1, true, 3.0786, 1.9242, 0, "Paine & Hann" }},                    
    { 316, { "red maple", "Acer", 2, false, 2.17, 0.491, 0, "Russell & Weiskittel " }},                 
    { 351, { "red alder", "Alnus", 1, true, 8, 1.53, 0, "Smith" }},                                     
    { 361, { "Pacific madrone", "Arbutus", 1, true, 3.4299, 1.3532, 0, "Paine & Hann" }},               
    { 371, { "yellow birch", "Betula", 2, false, 4.04, 0.308, 0, "Russell & Weiskittel " }},            
    { 375, { "paper birch", "Betula", 2, false, 1.48, 0.623, 0, "Russell & Weiskittel " }},             
    { 379, { "gray birch", "Betula", 2, false, 2.24, 0.382, 0, "Russell & Weiskittel " }},              
    { 431, { "golden chinkapin", "Castanopsis", 1, true, 2.9794, 1.5512, -0.014161, " Paine & Hann" }}, 
    { 531, { "American beech", "Fagus", 2, false, 2.93, 0.434, 0, "Russell & Weiskittel " }},           
    { 631, { "tanoak", "Lithocarpus", 1, true, 4.4443, 1.704, 0, "Paine & Hann" }},                     
    { 746, { "quaking aspen", "Populus", 2, false, 1.31, 0.586, 0, "Russell & Weiskittel " }},          
    { 747, { "black cottonwood", "Populus", 1, true, 0.5, 1.62, 0, "Smith" }},                          
    { 768, { "bitter cherry", "Prunus", 1, true, 8, 1.53, 0, "Use RA (Smith)" }},                       
    { 805, { "canyon live oak", "Quercus", 1, true, 5, 1.69, 0, "Warbington & Levitan" }},              
    { 815, { "Oregon white oak", "Quercus", 1, true, 3.0786, 1.9242, 0, "Paine & Hann" }},              
    { 818, { "California black oak", "Quercus", 1, true, 3.3625, 2.0303, -0.0073307, " Paine & Hann" }},
    { 833, { "northern red oak", "Quercus", 2, false, 4.08, 0.31, 0, "Russell & Weiskittel " }},        
    { 920, { "willow", "Salix", 1, true, 8, 1.53, 0, "Use RA (Smith)" }},   
}};


std::vector<double> compute_mcw( const std::vector<int> &fia,
                                 const std::vector<double> &dbh,
                                 const bool imperial_units,
                                 const int default_fia )
{
    MCWPARMS mcw_data;
    size_t n = fia.size();
    std::vector<double> mcw_hat( n, 0.0 );

    for( size_t i = 0; i < n; ++i )
    {
        try {
            mcw_data = mcw_parms.at( fia[i] );
        } catch( std::out_of_range &e ) {
            try {
                mcw_data = mcw_parms.at( default_fia );
            } catch( std::out_of_range &e2 ) {
                throw( "No MCW data for supplied species and default species.\n" );
            }
        }

        double d = dbh[i];
        if( imperial_units && !mcw_data.imperial_units ) {
            d *= 2.54;
        } else if( !imperial_units && mcw_data.imperial_units ) {
            d /= 2.54;
        }

        mcw_hat[i] = mcw_data.EQ_type == 1 ?  
                        mcw_data.A + (mcw_data.B * dbh[i]) + (mcw_data.C * dbh[i] * dbh[i]) : 
                        mcw_data.A * std::pow(dbh[i], mcw_data.B);

        if( imperial_units && !mcw_data.imperial_units ) {
            mcw_hat[i] *= 3.2808;
        } else if( !imperial_units && mcw_data.imperial_units ) {
            mcw_hat[i] /= 3.2808;
        }
    }

    return mcw_hat;

}

std::vector<std::pair<int, std::string>> get_mcw_species()
{
    std::vector<std::pair<int, std::string>> species;
    for( auto &[fia, m] : mcw_parms )
    {
        species.push_back( std::pair( fia, m.name ) );
    }

    return species;
}