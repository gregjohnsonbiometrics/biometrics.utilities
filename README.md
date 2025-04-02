# biometrics.utilities
A computationally efficient implementation in R[^1] of tree and plot level metrics used in Forest Biometrics.

## Purpose
This R package is a collection of metrics (mostly focused on competition indicies) that are useful in describing
tree and stand conditions. Several of them are somewhat complex and computationally expensive. The code here
translates these metrics into C++ for computational efficiency and hopefully sound implementations of the 
concepts.

This is a collaborative effort with David Marshall.

If you have suggestions for additional metrics, let me know ([Greg Johnson](mailto:greg@nosnhoj.org)).

## Metrics

### Distance Independent

* Basal Area ([ba]())
* Crown Closure at Tree Tip ([cch](#crown-closure-at-tree-tip-cch))
* Basal Area in Larger Trees ([bal](#basal-area-in-larger-trees-bal))
* Crown Competition Factor in Larger Trees ([ccfl](#crown-competition-factor-in-larger-trees-ccfl))
* Dominant Height ([dominant_height](#dominant-height))
* Quadratic Mean Diameter ([qmd](#quadratic-mean-diameter-qmd))
* Wilson's Relative Spacing ([relative_spacing](#wilsons-relative-spacing-relative_spacing))
* Crown Competition Factor ([ccf](#crown-competition-factor-ccf))
* Curtis' Relative Density ([curtis_rd](#curtis-relative-density-curtis_rd))
* Reineke's Stand Density Index ([reineke_sdi](#reinekes-stand-density-index-reineke_sdi))
* Glover and Hool Index ([Glover_Hool](#glover-and-hool-competition-index-glover_hool))


### Distance Dependent

* Clark-Evans R Aggregation Index ([Clark_Evans_R](#clark-evans-aggregation-index-clark_evans_r))
* Hegyi's Distance-weighted size ratio ([Hegyi](#hegyi-competition-index-hegyi))
* Arney's Competitive Stress Index ([Arney_CSI](#arneys-competitive-stress-index-arney_csi))
* Area Potentially Available ([APA](#area-potentially-available-apa))

### Utilities

* Maximum Crown Width Estimate ([mcw](#maximum-crown-width-mcw))
* Fit Height-DBH Curves ([hd_fit](#fit-height---dbh-curves-hd_fit))
* Predict Height from DBH ([hd_predict](#fit-height---dbh-curves-hd_fit))

-------------

### Basal Area (`ba`)

`ba( dbh, expansion, imperial_units )`

where:
- `dbh` = diameter at breast height (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`ba()` is a convenience function to compute basal area per unit area (square feet per acre or square meters per hectare depending
on the value of `imperial_units`).

$ba = \sum( dbh_i^2 \times expansion_i \times k)$

where k converts squared diameters to square feet per acre or meters per hectare depending on `imperial_units`.

`ba` returns a double containing the basal area.

-------------

### Crown closure at tree tip (`cch`)

`cch( species, dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units )`

where:
- `species` = species code (vector)
- `dbh` = diameter at breast height (vector)
- `height` = total height (vector)
- `dacb` = distance above crown base (vector)
- `lcw` = largest crown width (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `parameters` = `data.frame` of 3 parameters for each species used in the `cwa` (crown width above) equation
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`cch`[^2] has been used in growth models to quantify crown competition (typically in conifers). `cch` is the crown area as a fraction of an acre or hectare (depending on units of measure) that lies in a horizontal plane tangential to the tip of the tree. A value of 1.0 for a tree would indicate that the other trees in the stand or plot have crown area at the tip of the subject tree that completely covers the area. In contrast to other existing implementations, this function *does not* use an interpolation table scheme.

`dacb` estimates the distance between the crown base to the largest crown width point in the tree. Hann (1999)[^3] found that this distance was proportional to crown length (`cl`).

`lcw` is an estimate for each tree of the largest crown width. This is usually estimated from species-specific prediction equations[^3].

The `cwa` equation is used to compute the crown width at a point above the largest crown width. The point is usually described by its relative position (i.e. 0 = tip of the tree, 1 = at the largest crown width point). The form of the equation is:

$cwa =  rp^{(\beta_0 + \beta_1 rp^{0.5} + \beta_2 height / dbh)}$

where: rp = relative position, and $\beta_0$ - $\beta_2$ are parameters to be supplied by the user in the `parameters` `data.frame`. The `data.frame` has the
following members:

- species : species code
- b0 : $\beta_0$ parameter
- b1 : $\beta_1$ parameter
- b2 : $\beta_2$ parameter

If a species is missing from the `data.frame`, a cone is used (1.0, 0.0, 0.0).

`cch` returns a vector of crown closures for each tree.

-------------

### Basal Area in Larger Trees (`bal`)

`bal( dbh, expansion, imperial_units)`

where:

- `dbh` = diameter at breast height (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

The `bal` function sorts in-place the dbh and expansion vectors in decreasing order of dbh and computes cumulative basal area in larger trees. Ties in dbh are handled appropriately (each tree has the same `bal`). Ties in `dbh` are handled appropriately (each tree has the same `bal`). Thus the largest tree (by `dbh`) has a `bal` of 0.0 and the smallest tree has a `bal` equal to the stand or plot basal area less the basal area of the subject tree.

The `bal` function returns a vector in the original tree order of the input variables.

`bal` is typically used to quantify one-sided competition in diameter growth equations.

-------------

### Crown Competition Factor in Larger Trees (`ccfl`)

`ccfl( dbh, expansion, imperial_units)`

where:

- `dbh` = diameter at breast height (vector)
- `mcw` = maximum crown width of open grown tree (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

The `ccfl` function computes the [crown competition factor](#crown-competition-factor-ccf) for a tree based on a subset of trees in the stand or plot that have a `dbh` greater than the subject tree. The function sorts the `dbh`, `mcw`, and `expansion` vectors in-place in decreasing order of `dbh` and computes cumulative maximum crown area in larger trees. Ties in `dbh` are handled appropriately (each tree has the same `ccfl`).

`mcw` is an estimate of the maximum crown width (open grown) for a tree of diameter `dbh`. Hann (1999)[^3] used the following equation form:

$mcw = \beta_0 + \beta_1 dbh + \beta_2 dbh^2$

for each species in Southwest Oregon.

The `ccfl` function returns a vector in the original tree order of the input variables.

`ccfl` is often used to quantify one-sided competition in growth equations.

-------------

### Dominant Height

`dominant_height( height, dbh, expansion, dominant_cohort_size, method )`

where:

- `height` = total height (vector)
- `dbh` = diameter at breast height (vector)
- `expansion` = expansion factors (vector) 
- `dominant_cohort_size` = number of trees in the dominant height cohort
- `method` = 0 (default), 1, or 2 (see below for definitions)

`dominant_height` computes the expansion factor weighted average height of trees in the defined dominant tree cohort. Each method defines the cohort differently:

- 0 = average height of the `dominant_cohort_size` trees by decreasing `dbh`
- 1 = average height of the `dominant_cohort_size` trees by decreasing `height`
- 2 = Lorey height (height of the tree of average basal area)

Common values for `dominant_cohort_size` are 100 and 40 trees. The cohort size is ignored and should be 0.0 for Lorey height (option 2).

`dominant_cohort_size` returns a scalar with the specified dominant height for the stand or plot.

-------------

### Quadratic Mean Diameter (qmd)

`qmd( dbh, expansion )`

where:

- `dbh` = diameter at breast height (vector)
- `expansion` = expansion factors (vector) 

`qmd` computes the quadratic mean diameter[^8] (basal area weighted mean diameter) for the stand or plot.

`qmd` returns a scalar with the quadratic mean diameter.

-------------

### Wilson's Relative Spacing (`relative_spacing`)

`relative_spacing( expansion, dom_height, imperial_units )`

where:

- `expansion` = expansion factors (vector)
- `dom_height` = dominant height (scalar) 
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`relative_spacing` computes Wilson's Relative Spacing[^7]. It is defined as the average tree spacing relative to dominant height. Smaller
values indicate more crowding or inter-tree competition.

`relative_spacing` returns a scalar with the spacing as a fraction of dominant height.

-------------

### Crown Competition Factor (`ccf`)

`ccf( crown_width, expansion, imperial_units )`

where:

- `crown_width` = open-grown crown widths (vector) 
- `expansion` = vexpansion factors (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`ccf` computes the Crown Competition Factor[^4] for the stand or plot. `ccf` is the ratio of the open-grown crown area of all trees as a percentage of an acre (or hectare depending on `imperial_units`). In the imperial units case:

$ccf = 100 \frac{\sum{(ca_i expansion)}}{43560}$

where $ca_i$ is the open-grown crown area for tree $i$.

`ccf` returns a scalar with the crown competition factor for the stand or plot.

-------------

### Curtis' Relative Density (`curtis_rd`)

`curtis_rd( dbh, expansion, imperial_units )`

where:

- `dbh` = diameter at breast height (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`curtis_rd` compute Curtis' Relative Density[^5].

Curtis' relative density measure is computed as:

$rd = G/(Dg^{0.5})$

where G is basal area and Dg is quadratic mean stand diameter.

`curtis_rd` returns a scalar with Curtis' Relative Density.

-------------

### Reineke's Stand Density Index (`reineke_sdi`)

`reineke_sdi( dbh, expansion, imperial_units )`

where:

- `dbh` = diameter at breast height (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`reineke_sdi` computes Reineke's Stand Density Index[^6]. Reineke developed a stand density index (SDI) that relates the current stand density to an equivalent density in a stand with a quadratic mean diameter (Dq) of 10 inches. SDI can be expressed for the imperial units case as:

$sdi = N(\frac{Dq}{10})^b$

where `sdi`= Reineke’s stand density index, N = trees per acre, Dq = quadratic mean diameter (inches), b = exponent of Reineke’s equation, often reported to equal –1.605.

`reineke_sdi` returns a scalar with the Reineke's Stand Density Index.

-------------

### Glover and Hool Competition Index (`Glover_Hool`)

`Glover_Hool( dbh, expansion, use_arithmetic, imperial_units )`

where:

- `dbh` = diameter at breast height (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `use_arithmetic` = boolean where TRUE = use arithmetic mean diameter, FALSE = use quadratic mean diameter
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`Glover_Hool` computes the Glover and Hool (1979)[^16] competition index. The index is interpreted as the ratio
of a tree's basal area to the basal area of the tree of mean diameter. Glover and Hool used
the arithmetic mean and a common variation is to use the quadratic mean (use the `use_arithmetic` flag to select 
the desired method).

The index $G_i$ is:

$G_i = dbh_i^2 / \overline{dbh}^2$

where: $dbh_i$ is the diameter of tree `i` and $\overline{dbh}$ is the mean diameter (either arithmetic or quadratic).


`Glover_Hool` returns a vector of competition indicies for each tree in their original order.

-------------

### Clark-Evans Aggregation Index (`Clark_Evans_R`)

`Clark_Evans_R( x, y, plotarea, plot_x, plot_y )`

`Clark_Evans_R_circle( x, y, plot_area, plot_center_x, plot_center_y, plot_radius )`

where:

- `x` = x coordinates of trees on a plot (vector)
- `y` = y coordinates of trees on a plot (vector)
- `plotarea` = area of the plot polygon or circle to use if polygon or circle coordinates are not supplied
- `plot_x` = vector of plot vertices x coordinates
- `plot_y` = vector of plot vertices y coordinates
- `plot_center_x` = x coordinate of plot center
- `plot_center_y` = y coordinate of plot center
- `plot_radius` = radius of plot

Clark_Cvans_R` computes the Clark and Evans (1954)[^9] aggregation index. The aggregation index R is a measure of clustering or ordering of trees on a plot. It is the ratio of the observed mean nearest neighbor distance in the trees to that expected for a Poisson point process of the same intensity. A value R > 1 suggests ordering, while R < 1 suggests clustering (unequal inter-tree competition). R has been proposed as a two-sided, distance-dependent tree competition metric.

This implementation can use Donnelly's[^11] edge correction if polygon or circle coordinates are supplied; otherwise an uncorrected estimate is returned.

$R = \frac{\frac{\sum{ d_i }}{N}}{(\frac{A}{N})^{0.5}/2}$

`Clark_Evans_R` and `Clark_Evans_R_circle` return a scalar with the Clark Evans R statistic.

-------------

### Hegyi Competition Index (`Hegyi`)

`hegyi( x, y, dbh, plot_x, plot_y, imperial_units )`

where:

- `x` = x coordinates of trees on a plot (vector)
- `y` = y coordinates of trees on a plot (vector)
- `dbh` = diameter at breast height (vector)
- `plot_x` = vector of plot vertices x coordinates
- `plot_y` = vector of plot vertices y coordinates
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`hegyi` computes the distance-weighted size ratio (a two-sided competition index) based on Hegyi (1974)[^10] for trees within a 6-meter fixed radius plot.

The ratio for the ith tree in the 6-meter radius plot is:

$hegyi_i = \sum{\frac{ba_j/ba_i}{d_{ij}}}$

where $ba_j$ and $ba_i$ are the basal areas of the jth and ith tree respectively, $d_{ij}$ is the distance between tree i and j.

Trees from a plot of arbitrary size can be used. The Hegyi ratio for each tree will be computed based on its neighbors within the 6-meter boundary. If `imperial_units`
is TRUE, the coordinates will be converted to meters prior to calculations.

This version adjusts for edge effects using Ripley's[^14] edge correction if plot coordinates are supplied. If no plot corners are available, the index is unadjusted.

`hegyi` returns ratios for each tree in the input vectors (preserving their order).

-------------

### Arney's Competitive Stress Index (`Arney_CSI`)

`Arney_CSI( x, y, dbh, mcw )`

where:

- `x` = x coordinates of trees on a plot (vector)
- `y` = y coordinates of trees on a plot (vector)
- `dbh` = diameter at breast height (vector)
- `mcw` = maximum crown width of open grown tree (vector)

`Arney_CSI` computes Arney's (1973)[^12] Competitve Stress Index (`CSI`). CSI is the sum of the percentage competing trees crown area overlaping
a subject tree to the subject tree's crown area. `CSI` is essentially a transformation of Gerrard's (1969)[^13] Competition Quotient using maximum crown area for 
each tree as its competition circle (Gerrard used an empirically derived circle radius dependent on `dbh`).

$CSI_i = 100 \sum{\frac{AO_j}{CA_i}}$

where $AO_j$ is the area of overlap of tree j on subject tree i, $CA_i$ is the crown area of the subject tree i, and $CSI_i$ is the
competitive stress index for tree i.

This version currently does not adjust for edge effects.

`Arney_CSI` returns a vector of CSI values for each tree in the original order.

-------------

### Maximum Crown Width (`mcw`)

`mcw( fia, dbh, imperial_units, default_fia )`

where:

- `fia` = FIA numeric species code (vector)
- `dbh` = diameter at breast height (vector)
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric
- `default_fia` = FIA species code to use if supplied `fia` code is not in parameter table.

`mcw` uses publicly available parameter estimates for one of two maximum crown width estimatation functions. The equations are either two or three parameter
equations of the following forms:

1: $\widehat{mcw} = A + B dbh + C dbh^2$, or

2: $\widehat{mcw} = A dbh^B$

`mcw` returns the maximum crown width for each tree in the input vectors in feet or meters depending on `imperial_units`.

A helper function: `mcw_species()` is available to produce a `data.frame` of FIA species codes and species names with parameter estimates available.

-------------

### Fit Height - DBH Curves (`hd_fit`)

`hd_fit( fia, dbh, height, bh )`

`hd_predict( fit, fia, dbh, bh )`

where:

- `fia`    = FIA numeric species code (vector)
- `dbh`    = diameter at breast height (vector)
- `height` = total height (vector)
- `bh`     = height to breast height (scalar)
- `fit`    = a height-dbh model fit returned by the `hd_fit` function.

`hd_fit` is a function to fit height-dbh curves of one of the forms found in Curtis (1967)[^15]:

$height = bh + e^{(\beta_0 + \beta_1 dbh^{\beta_2})}$

where $\beta$ s are parameters to be estimated.

Separate curves are fit to each species. If a species has less than 3 observations, a parameter set of 0.0 is generated.

Fitting height-dbh curves is often difficult due to limited measurement data (either in observation count or in range, or both). We are using David Marshall's 
technique of fitting a linearized form of the equation while iterating over a range of $\beta_2$ values (-0.1 to -1.0). The $\beta_2$ value yielding the lowest sum of squared errors
(SSE) is chosen.

`hd_fit` returns a `data.frame` of parameter estimates by FIA species code.

`hd_predict` is used to predict `height` for a vector of `dbh` values given parameter estimates from `hd_fit`.

-------------

### Area Potentially Available (`APA`)

`APA( x, y, dbh, poly_x, poly_y, weighted )`

`APA_Polygons( tree_id, x, y, dbh, poly_x, poly_y, weighted )`

where:

- `tree_id` = unique identification number for each tree (vector)
- `x` = x coordinates of trees on a plot (vector)
- `y` = y coordinates of trees on a plot (vector)
- `dbh` = diameter at breast height (vector)
- `poly_x` = x coordinates of the plot bottom left and top right corners
- `poly_y` = y coordinates of the plot bottom left and top right corners
- `weighted` = boolean flag to select unweighted (`false`) for a standard Voronoi polygon area, or weighted (`true`) for a Voronoi polygon constructed by using `dbh` as a weight.

`APA` computes the Area Potentially Available (APA) using Brown's (1965)[^17] method. Polygons are constructed around the subject tree
by the intersection of the perpendicular bisectors of the distance between the subject tree and competitors (creating a 
Voronoi tesselation of the plot).

`APA_Polygons` builds the APA polygons and returns a `data.frame` with the polygon coordinates.

This version currently does not adjust for edge effects.

This version currently does not implement the `weighted` option.

The Voronoi tesselation is computed using Fortune's algorithm[^18] from code derived from Pierre Vigier Copyright (C) 2018, and is provided under GNU Lesser General Public License (see <http://www.gnu.org/licenses/>).

`APA` returns a vector of APA values (in square units of measure used for the coordinates) for each tree in their original order.

## R Packages

The source and binary packages can be found in the repository:

- Windows Binary: [biometrics.utilities_1.0.zip](./biometrics.utilities_1.0.zip)
- Source: Windows Binary: [biometrics.utilities_1.0.tar.gz](./biometrics.utilities_1.0.tar.gz)

[^1]: R Core Team (2024). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

[^2]: Ritchie, M.W. and D.W. Hann. 1990. Equations for predicting height growth of six conifer species in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Paper 54. 12p.

[^3]: Hann, D.W. 1999. An adjustable predictor of crown profile for stand-grown Douglas-fir trees. For. Sci. 45: 217–225.

[^4]: Krajicek, J.E., K.A. Brinkman, and F.S. Gingrich.  1961.  Crown competition: a measure of density.  For. Sci. 7:36 – 42. 

[^5]: Curtis, R.O. 1982. A Simple Index of Stand Density for Douglas-fir. Forest Sci., Vol. 28, No. 1, pp. 92-94.

[^6]: Reineke, L.H. 1933. Perfecting a stand density index for even-aged forests. Jour. Agric. Res.  46: 627 – 638. 

[^7]: Wilson, F.G. 1946. Numerical expression of stocking in terms of height. Journal of Forestry, 44:758–761.  

[^8]: Curtis, Robert O.; Marshall, David D. 2000. Why quadratic mean diameter? Western Journal of Applied Forestry, 15 (3): 137–139.

[^9]: Clark, P. J., & Evans, F. C. 1954. Distance to Nearest Neighbor as a Measure of Spatial Relationships in Populations. Ecology, 35(4), 445–453. <https://doi.org/10.2307/1931034>

[^10]: Hegyi, F. 1974. A simulation model for managing jack-pine stands. J. Fries (Ed.), Growth Models for Tree and Stand Simulation, Royal College of Forestry, Stockholm, Sweden (1974), pp. 74-90.

[^11]: Donnelly, K. 1978. Simulations to determine the variance and edge-effect of total nearest neighbour distance. In I. Hodder (ed.) Simulation studies in archaeology, Cambridge/New York: Cambridge University Press, pp 91–95.

[^12]: Arney, J.D. 1973. Tables for quantifying competitive stress on individual trees. Pacific Forest Research Centre, Canadian Forest Service, Victoria, BC. Information Report BC-X-78. 47p.

[^13]: Gerrard, D. J. 1969. Competition Quotient: a new measure of the competition affecting individual forest trees. Michigan State University, Agr. Exp. Sta. Res. Bull. 20. 32pp.

[^14]: Ripley, B.D. (1977) Modelling spatial patterns (with discussion). Journal of the Royal Statistical Society, Series B, 39, 172 -- 212.

[^15]: Curtis, R. 0. 1967. Height-diameter, and height-diameter-age equations for second growth Douglas-fir. For. Sci. 365-375.

[^16]: Glover G.R. and Hool, J.N. 1979. A basal area ratio predictor of loblolly pine plantation mortality. For. Sci. 25:275-282.

[^17]: Brown, G.S. 1965. Point density in stems per acre. N.Z. Forest Res. Notes No. 38.

[^18]: Fortune, S. (1987). A sweepline algorithm for Voronoi diagrams. Algorithmica, 2(1), 153–174. https://doi.org/10.1007/BF01840357