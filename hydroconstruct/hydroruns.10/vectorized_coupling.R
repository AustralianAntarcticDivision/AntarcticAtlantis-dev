library(angstroms)
library(rbgm)
library(dplyr)
delif <- function(x) {
  f <- filename(x)
  if (file.exists(f)) {
    ff <- gsub("grd$", "gri", f)
    if (file.exists(ff)) file.remove(ff)
  }
  invisible(NULL)
}

face_area <- function(bgm) {
  
}
## TODO
## DONE rotate u/v vectors on import (see ROMS details)
## check hyper diffusion scaling - I think this can be ignored
## check vertical flux calcs - I think this is ok

wrap180 <- function(x, lmin = -180) (x - lmin) %% 360 + lmin
sgn1 <- function(x) ifelse(x < 0, -1, 1)


cpolar <- raadtools:::cpolarfiles() %>% dplyr::filter(stringr::str_detect(fullname, "his_31"))
file_db <- purrr::map_df(cpolar$fullname, 
                         function(x) tibble::tibble(fullname = x,
                                                    time = tidync::tidync(x)$transforms$ocean_time$ocean_time) %>% 
                           dplyr::mutate(dim4_slice = row_number()))

## GRRR 
#time_steps <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + file_db$time
## ignore the file time steps and assume they daily from the start value
file_db$time_steps <-  seq(ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + file_db$time[1], by = "1 day", length.out = nrow(file_db))
time_step <- 86400

## get a BGM and read it
bfile <- bgmfiles::bgmfiles("antarctica_28")
bgm <- rbgm:::read_bgm(bfile)

roms_file <- file_db$fullname[1]
roms_ll <- romscoords(roms_file, transpose = TRUE)
boxes <- boxSpatial(bgm)
## width and height of each box
e <- vector("list", nrow(boxes))
for (i in seq_along(e)) {ex <- raster::extent(boxes[i, ]); e[[i]] <- c(xmin = xmin(ex), xmax = xmax(ex), ymin = ymin(ex), ymax = ymax(ex))}
boxes$width <- apply(do.call(rbind, e)[,1:2], 1, diff)
boxes$height <- apply(do.call(rbind, e)[,3:4], 1, diff)
faces <- faceSpatial(bgm)
faces_ll <- spTransform(faces, "+init=epsg:4326")
boxes_ll <- spTransform(boxes, "+init=epsg:4326")

faces$.x0 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][1, 1])
faces$.x1 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][2, 1])
faces$.y0 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][1, 2])
faces$.y1 <- purrr::map_dbl(coordinates(faces_ll), ~.x[[1]][2, 2])

romsbox <- romsmap(boxes, roms_ll)
roms_ll <- crop(roms_ll , romsbox, snap = "out")
romsface <- romsmap(faces, romscoords(roms_file, transpose = TRUE))

## dummy grids
u <- crop(romsdata3d(roms_file, varname = "u", slice = 1), romsbox)
v <- crop(romsdata3d(roms_file, varname = "v", slice = 1), romsbox)
temp <- crop(romsdata3d(roms_file, varname = "temp", slice = 1), romsbox)
salt <- crop(romsdata3d(roms_file, varname = "salt", slice = 1), romsbox)
vert <- subset(crop(romsdata3d(roms_file, varname = "w", slice = 1), romsbox), 1:31)
h <- crop(romsdata(roms_file, varname = "h"), romsbox, snap = "out")
hhh <- crop(romsdepth(roms_file), romsbox, snap = "out")
## snap h onto the first layer (the bottom) and calculate delta_rho (dz)
delta_rho <- calc(brick(-h, hhh), fun = function(x) diff(x))



deepest_depth <- max(raster::extract(h, tabularaster::cellnumbers(h, sf::st_as_sf(romsbox))$cell_), na.rm = TRUE)
## here we must use the same as we used for the nc output files
dz <- rev(rbgm::build_dz(-deepest_depth, zlayers = c(-Inf,  -2000, -1000, -750, -400, -300, 
                                                     -200, -100, -50, -20, 0)))
atlantis_depths <- head(-cumsum(c(0, dz)), -1)


rvalues <- function(x) {
  if (nchar(raster::filename(x)) > 0) x <- readAll(x)
  as.vector(raster::values(x))
}
romstab <- tabularaster::as_tibble(u, dim = FALSE) %>% 
  dplyr::rename(u = cellvalue) %>% dplyr::mutate(v = rvalues(v), 
                                                 temp = rvalues(temp), 
                                                 salt = rvalues(salt), 
                                                 vert = rvalues(vert),
                                                 depth = rvalues(hhh), 
                                                 #dz = rvalues(delta_rho),  ## not this, this is ROMS cell depth, not face delta
                                                # harea = rvalues(harea),
                                                 layer = findInterval(-depth, -atlantis_depths), 
                                                delta_layer = dz[layer])
delif(delta_rho)
#delif(harea)


box_index <- tabularaster::cellnumbers(u, romsbox)
box_index$boxid <- romsbox$box_id[box_index$object_]
## keep botz as we use it to filter too-low ROMS data
box_index$maxlayer  <- findInterval(-romsbox$botz[box_index$object_], -atlantis_depths)
face_index <- tabularaster::cellnumbers(u, romsface)
face_index$faceid <- romsface$.fx0[face_index$object_]
face_index$object_ <- NULL
face_index <- inner_join(face_index, faces@data[c("length", ".fx0", ".x0", ".y0", ".x1", ".y1")], c("faceid" = ".fx0"))
face_index <- face_index %>% 
  mutate(angle_face = (180 * atan2(.x0 - .x1, .y0 - .y1))/pi)   # angle in radians of tan(y/x) / pi resulting in angle in degrees
face_index$.x0 <- NULL
face_index$.x1 <- NULL
face_index$.y0 <- NULL
face_index$.y1 <- NULL


### FOR TRANSPORT NC FILE
transp_filename <- sprintf("%s_transport.nc", tempfile())
mass_filename  <- sprintf("%s_mass.nc", tempfile())
transp_params <- ""
mass_params <- ""

angle <- crop(angstroms:::romsangle(roms_file), romsbox)




angstroms:::create_transport(transp_filename, model_title = "Transport file Antarctica_28", bgmfilepath = bfile, 
                             bgmlevels = atlantis_depths, time_steps = file_db$time_steps)

angstroms:::create_mass(mass_filename, model_title = "Mass file Antarctica_28", bgmfilepath = bfile, 
                        bgmlevels = atlantis_depths, time_steps = file_db$time_steps)


vertical <- salinity <- temperature <- angstroms::rawdata(mass_filename, "temperature")
transport <- angstroms::rawdata(transp_filename, "transport")

# system.time(doitfun(1))
# user  system elapsed 
# 11.264   2.105  13.465 
# system.time(doitfun(300))
# user  system elapsed 
# 10.214   1.342  16.637 

doitfun <- function(itime) {
  roms_file <- file_db$fullname[itime]
  roms_slice <- file_db$dim4_slice[itime]
  u0 <- crop(romsdata3d(roms_file, varname = "u", slice = roms_slice), romsbox)
  v0 <- crop(romsdata3d(roms_file, varname = "v", slice = roms_slice), romsbox)
  
  ## apply required ROMs rotation  https://www.myroms.org/forum/viewtopic.php?f=3&t=295
  ul <- vector("list", nlayers(u0))
  vl <- vector("list", nlayers(u0))
  for (j in seq_along(ul)) {
    uv <- romsrotate(brick(u0[[j]], v0[[j]]), angle)
    ul[[j]] <- uv[[1]]
    vl[[j]] <- uv[[2]]
  }
  u0 <- brick(ul)
  v0 <- brick(vl)
  for (j in seq_along(ul)) {
    delif(ul[[j]])
    delif(vl[[j]])
  }
  rm(ul, vl)
  temp0 <- crop(romsdata3d(roms_file, varname = "temp", slice = roms_slice), romsbox)
  salt0 <- crop(romsdata3d(roms_file, varname = "salt", slice = roms_slice), romsbox)
  vert0 <- subset(crop(romsdata3d(roms_file, varname = "w", slice = roms_slice), romsbox), 1:31)
  
  romstab$u <- rvalues(u0)
  romstab$v <- rvalues(v0)
  romstab$temp <- rvalues(temp0)
  romstab$salt <- rvalues(salt0)
  romstab$vert <- rvalues(vert0)
  delif(u0)
  delif(v0)
  delif(temp0)
  delif(salt0)
  delif(vert0)
  
box_summ <- romstab %>% 
  inner_join(box_index[c("cell_", "boxid", "maxlayer")], c("cellindex" = "cell_")) %>% 
  dplyr::mutate(layer = pmin(layer, maxlayer)) %>%  ## push lower data into lowest layer
  group_by(boxid, layer) %>% summarize(temp = mean(temp, na.rm = TRUE), 
                                       salt = mean(salt, na.rm = TRUE), 
                                       vert = sum(vert, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(boxid, layer)


face_summ <- romstab %>% dplyr::select(-temp, -salt) %>% 
  inner_join(face_index, c("cellindex" = "cell_")) %>% 
  inner_join(box_index[c("cell_", "maxlayer")], c("cellindex" = "cell_")) %>% 
  dplyr::mutate(layer = pmin(layer, maxlayer)) %>%  ## push lower data into lowest layer
 ## mutate(angle_face = (180 * atan2(.x0 - .x1, .y0 - .y1))/pi)  %>%   # angle in radians of tan(y/x) / pi resulting in angle in degrees
  mutate(angle_flow = (180 * atan2(u, v))/pi) %>%
  mutate(relative_angle = wrap180(angle_flow - angle_face), 
         velocity_magnitude = sqrt(u^2 + v^2) * sgn1(relative_angle)) %>% 
  group_by(faceid, layer) %>% 
  ## mean because (i)  https://github.com/AustralianAntarcticDivision/EA-Atlantis-dev/issues/5#issuecomment-457453852
  summarize(flow = mean(velocity_magnitude * length * dz[layer], na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(faceid, layer)

 list(face = face_summ[c("faceid", "layer", "flow")], 
      box = box_summ[c("boxid", "layer", "temp", "salt", "vert")])
}


library(future)
plan(multicore)  ## use multiprocess on Windows
x <- future.apply::future_lapply(seq_len(nrow(file_db)), doitfun)


for (itime in seq_along(x)) {
  face_summ <- x[[itime]]$face
  box_summ <- x[[itime]]$box
  ## FIXME loop over these for itime
  fidx <- cbind(face_summ$layer, face_summ$faceid + 1, itime)
  bidx <- cbind( box_summ$layer, box_summ$boxid + 1, itime)
  temperature[bidx] <- box_summ$temp
  salinity[bidx] <- box_summ$salt
  vertical[bidx] <- box_summ$vert
  transport[fidx] <-face_summ$flow
}





library(ncdf4)


ncmass <- ncdf4::nc_open(mass_filename, write = TRUE)
nctran <- ncdf4::nc_open(transp_filename, write = TRUE)


ncvar_put(nctran, "transport", transport)
ncvar_put(ncmass, "temperature", temperature)
ncvar_put(ncmass, "salinity", salinity)
ncvar_put(ncmass, "verticalflux", vertical)
nc_close(nctran)
nc_close(ncmass)

file.copy(mass_filename, ".")
file.copy(transp_filename, ".")
