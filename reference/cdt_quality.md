# Triangle quality measures

Per-triangle area, angles, edge lengths and radius-edge ratio, for
judging a mesh or checking that refinement met its targets.

## Usage

``` r
cdt_quality(x, T = NULL)
```

## Arguments

- x:

  a result of
  [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
  (a list with `P` and `T`), or a two-column matrix of vertex
  coordinates, in which case `T` is needed

- T:

  triangle vertex indices (m x 3, 1-based), when `x` is a matrix

## Value

a data frame with one row per triangle and columns `area` (always
non-negative), `min_angle` and `max_angle` (degrees), `min_edge` and
`max_edge`, and `radius_edge` (circumradius over shortest edge; the
Ruppert bound on it is about 1.41 for a 20.7 degree minimum angle, and
`Inf` for a degenerate triangle). Use
[`summary()`](https://rdrr.io/r/base/summary.html) on the result for the
distributions.

## Examples

``` r
## an irregular polygon with a sharp corner, refined by area and angle
x <- c(0, 2, 2, 1, 0)
y <- c(0, 0, 2, 0.4, 2)
r <- cdt_triangulate(x, y, 1:5, c(2:5, 1), max_area = 0.05, min_angle = 25)
q <- cdt_quality(r)
summary(q[c("area", "min_angle", "radius_edge")])
#>       area            min_angle      radius_edge    
#>  Min.   :0.008281   Min.   :26.48   Min.   :0.5979  
#>  1st Qu.:0.027521   1st Qu.:38.21   1st Qu.:0.7071  
#>  Median :0.031250   Median :42.81   Median :0.7358  
#>  Mean   :0.031579   Mean   :41.90   Mean   :0.7672  
#>  3rd Qu.:0.039806   3rd Qu.:45.00   3rd Qu.:0.8084  
#>  Max.   :0.049642   Max.   :56.75   Max.   :1.1214  
## both refinement targets were met
c(max_area = max(q$area), min_angle = min(q$min_angle))
#>    max_area   min_angle 
#>  0.04964213 26.47852647 
```
