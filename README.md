# biometrics.utilities
A computationally efficient implementation in R[^1] of tree and plot level metrics used in Forest Biometrics.

## Purpose
This R package is a collection on metrics (mostly focused on competition indicies) that are useful in describing
tree and stand conditions. Several of them are somewhat complex and computationally expensive. The code here
translates these metrics into C++ for computational efficiency and hopefully sound implementations of the 
concepts.

If you have suggestions for additional metrics, let me know ([Greg Johnson](mailto:greg@nosnhoj.org)).

## Metrics

* Crown Closure at Tree Tip ([cch](#crown-closure-at-tree-tip-cch))
* Basal Area in Larger Trees ([bal](#basal-area-in-larger-trees-bal))
* Crown Competition Factor in Larger Trees ([ccfl](#crown-competition-factor-in-larger-trees-ccfl))
* Dominant Height ([dominant_height](#dominant-height))
* Quadratic Mean Diameter ([qmd](#quadratic-mean-diameter-qmd))
* Wilson's Relative Spacing ([relative_spacing](#wilsons-relative-spacing-relative_spacing))
* Crown Competition Factor ([ccf](#crown-competition-factor-ccf))
* Curtis' Relative Density ([curtis_rd](#curtis-relative-density-curtis_rd))
* Reineke's Stand Density Index ([reineke_sdi](#reinekes-stand-density-index-reineke_sdi))

-------------

### Crown closure at tree tip (`cch`)

`cch( dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units )`

where:
- `dbh` = diameter at breast height (vector)
- `height` = total height (vector)
- `dacb` = distance above crown base (vector)
- `lcw` = largest crown width (vector)
- `expansion` = tree factor to expand to per area basis (vector)
- `parameters` = vector of 3 parameters used in the `cwa` (crown width above) equation
- `imperial_units` = boolean where TRUE is imperial, FALSE is metric

`cch`[^2] has been used in growth models to quantify crown competition (typically in conifers). `cch` is the crown area as a fraction of an acre or hectare (depending on units of measure) that lies in a horizontal plane tangential to the tip of the tree. A value of 1.0 for a tree would indicate that the other trees in the stand or plot have crown area at the tip of the subject tree that completely covers the area. In contrast to other existing implementations, this function *does not* use an interpolation table scheme.

`dacb` estimates the distance between the crown base to the largest crown width point in the tree. Hann (1999)[^3] found that this distance was proportional to crown length (`cl`).

`lcw` is an estimate for each tree of the largest crown width. This is usually estimated from species-specific prediction equations{^3].

The `cwa` equation is used to compute the crown width at a point above the largest crown width. The point is usually described by its relative position (i.e. 0 = tip of the tree, 1 = at the largest crown width point). The form of the equation is:

$cwa =  rp^{(\beta_0 + \beta_1 rp^{0.5} + \beta_2 height / dbh)}$

where: rp = relative position, and $\beta_0$ - $\beta_2$ are parameters to be supplied by the user in the `parameters` vector.

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

Sorts the `dbh`, `mcw`, and `expansion` vectors in-place in decreasing order of `dbh` and computes cumulative maximum crown area in larger trees. Ties in `dbh` are handled appropriately (each tree has the same `ccfl`).

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

Common values for `dominant_cohort_size` are 100 and 40 trees.

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
