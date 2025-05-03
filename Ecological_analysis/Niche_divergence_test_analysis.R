library(humboldt) # https://github.com/jasonleebrown/humboldt
library(raster)
library(dismo)
library(dplyr)

## 1. Load bioclimatic variables 
bio_tif = stack("./chelsa.v2.1.bio.tif")
names(bio_tif) = paste0("bio_", 1:19)
keep.bio = c("bio_10", "bio_4", "bio_11", "bio_17")
bio_tif = subset(bio_tif, keep.bio)

## 2. Load occurrence data for three species
pts_fratercula = read.csv("./occ_fratercula.filter.csv", header = TRUE)
pts_davidi     = read.csv("./occ_davidi.filter.csv", header = TRUE)
pts_hueti      = read.csv("./occ_hueti.filter.csv", header = TRUE)

## 3. Crop environmental data to study region and convert to data.frame
region.in <- extent(95, 122, 15, 34.5)
bio_tif_c = crop(bio_tif, region.in)
env.points <- rasterToPoints(bio_tif_c, fun = NULL, spatial = FALSE) %>% as.data.frame()

## 4. Format occurrence data for humboldt
pts_fratercula$sp = "fratercul"
pts_fratercul = pts_fratercula[, c("sp", "x", "y")]
pts_davidi$sp = "davidi"
pts_davidi = pts_davidi[, c("sp", "x", "y")]
pts_hueti$sp = "hueti"
pts_hueti = pts_hueti[, c("sp", "x", "y")]

### Pairwise comparisons using humboldt###

### fratercula vs hueti

# NOT: Full environmental space comparison (tests niche equivalency/divergence)
full_not <- humboldt.doitall(
  inname = "fratercula.hueti.not",
  env1 = env.points, env2 = env.points,
  sp1 = pts_fratercula, sp2 = pts_hueti,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 0, reductype = "PCA", non.analogous.environments = "YES",
  correct.env = TRUE, env.trim = TRUE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 100, b.reps = 100,
  nae = "YES", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = FALSE,
  run.silent = FALSE, ncores = 10
)

# NDT: Shared environmental space comparison (tests niche evolution/divergence)
shared_ndt <- humboldt.doitall(
  inname = "fratercula.hueti.ndt",
  env1 = env.points, env2 = env.points,
  sp1 = pts_fratercula, sp2 = pts_hueti,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO",
  correct.env = TRUE, env.trim = FALSE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 100, b.reps = 100,
  nae = "NO", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = FALSE,
  run.silent = FALSE, ncores = 10, color.ramp = 4
)

### fratercula vs davidi

# NOT: Full space
full_not_2 <- humboldt.doitall(
  inname = "fratercula.davidi.not",
  env1 = env.points, env2 = env.points,
  sp1 = pts_fratercula, sp2 = pts_davidi,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 0, reductype = "PCA", non.analogous.environments = "YES",
  correct.env = TRUE, env.trim = TRUE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 100, b.reps = 100,
  nae = "YES", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = FALSE,
  run.silent = FALSE, ncores = 10, color.ramp = 4
)

# NDT: Shared space
shared_ndt_2 <- humboldt.doitall(
  inname = "fratercula.davidi.ndt",
  env1 = env.points, env2 = env.points,
  sp1 = pts_fratercula, sp2 = pts_davidi,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO",
  correct.env = TRUE, env.trim = TRUE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 100, b.reps = 100,
  nae = "NO", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = TRUE,
  run.silent = FALSE, ncores = 8, color.ramp = 4
)

### hueti vs davidi

# NOT: Full space
full_not_3 <- humboldt.doitall(
  inname = "hueti.davidi.not",
  env1 = env.points, env2 = env.points,
  sp1 = pts_hueti, sp2 = pts_davidi,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 0, reductype = "PCA", non.analogous.environments = "YES",
  correct.env = TRUE, env.trim = TRUE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 100, b.reps = 100,
  nae = "YES", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = TRUE,
  run.silent = FALSE, ncores = 8, color.ramp = 4
)

# NDT: Shared space
shared_ndt_3 <- humboldt.doitall(
  inname = "hueti.davidi.ndt",
  env1 = env.points, env2 = env.points,
  sp1 = pts_hueti, sp2 = pts_davidi,
  rarefy.dist = 0, rarefy.units = "km", env.reso = 0.08333333,
  reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO",
  correct.env = TRUE, env.trim = TRUE, env.trim.type = "RADIUS",
  trim.buffer.sp1 = 50, trim.buffer.sp2 = 50,
  pcx = 1, pcy = 2, col.env = e.var, e.var = c(3:ncol(env.points)),
  R = 100, kern.smooth = "auto", e.reps = 300, b.reps = 300,
  nae = "NO", thresh.espace.z = 0.001,
  p.overlap = TRUE, p.boxplot = FALSE, p.scatter = TRUE,
  run.silent = FALSE, ncores = 8, color.ramp = 4
)
