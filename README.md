# biometrics.utilities
A computationally efficient implementation in R[^1] of tree and plot level metrics used in Forest Biometrics.

## Purpose
This R package is a collection on metrics (mostly focused on competition indicies) that are useful in describing
tree and stand conditions. Several of them are somewhat complex and computationally expensive. The code here
translates these metrics into C++ for computational efficiency and hopefully sound implementations of the 
concepts.

## Metrics

### Tree-level

#### Crown closure at tree tip (CCH)

cch( dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units )
where:
- dbh = diameter at breast height (vector)
- height = total height (vector)
- dacb = distance above crown base (vector)
- lcw = largest crown width (vector)
- expansion = tree factor to expand to per area basis (vector)
- parameters = vector of 3 parameters used in the `cwa` (crown width above) equation
- imperial_units = boolean where TRUE is imperial, FALSE is metric

CCH[^2] has been used in various growth models to quantify crown competition (typically in conifers). CCH is the crown area as a fraction of an acre or hectare (depending on units of measure) that lies in a horizontal plane tangential to the tip of the tree. A value of 1.0 for a tree would indicate that the other trees in the stand or plot have crown area at the tip of the subject tree that completely covers the area.

`dacb` estimates the distance between the crown base to the largest crown width point in the tree. Hann (1999)[^3] found that this distance was proportional to crown length (`cl`).

`lcw` is an estimate for each tree of the largest crown width. This is usually estimated from species-specific prediction equations{^3].

The `cwa` equation is used to compute the crown width at a point above the largest crown width. The point is usually described by its relative position (i.e. 0 = tip of the tree, 1 = at the largest crown width point). The form of the equation is:

$cwa =  rp^{(\beta_0 + \beta_1 rp^{0.5} + \beta_2 height / dbh)}$

where: rp = relative position, and $\beta_0$ - $\beta_2$ are parameters to be supplied by the user in the `parameters` vector.



[^1]: R Core Team (2024). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.

[^2]: Hann, D.W. and C. H. Wang. 1990. Mortality equations for individual trees in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Bulletin 67. 17p.

[^3]: [Hann, D.W. 1999. An adjustable predictor of crown profile for stand-grown Douglas-fir trees. For. Sci. 45: 217–225.]
