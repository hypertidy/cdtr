## Contours and a validation grid for the terrain vignette: kunanyi / Mount
## Wellington, Tasmania, from the Tasmania 2 m DEM (The LIST, (c) State of
## Tasmania, CC BY 3.0 AU), read as a COG through sds::tas_dem().
library(terra)

dsn <- sds::tas_dem()
r <- rast(dsn)
ctr <- c(519390, 5250707)  # GDA94 / MGA zone 55
e <- ext(ctr[1] - 3000, ctr[1] + 3000, ctr[2] - 3000, ctr[2] + 3000)
dem <- resample(crop(r, e), rast(e, resolution = 40, crs = crs(r)), method = "average")

levels <- seq(200, 1250, by = 50)
lines <- simplifyGeom(as.contour(dem, levels = levels), tolerance = 5)
g <- as.data.frame(geom(lines))
contours <- data.frame(line = as.integer(factor(paste(g$geom, g$part, sep = "."),
                                                levels = unique(paste(g$geom, g$part, sep = ".")))),
                       x = round(g$x, 1), y = round(g$y, 1),
                       z = lines$level[g$geom])

grid <- list(z = round(matrix(values(dem)[, 1], nrow(dem), ncol(dem), byrow = TRUE), 1),
             xmin = xmin(dem), xmax = xmax(dem), ymin = ymin(dem), ymax = ymax(dem),
             crs = "EPSG:28355")
saveRDS(list(contours = contours, grid = grid), "inst/extdata/tas_wellington.rds",
        compress = "xz")
