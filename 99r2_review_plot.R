# Review - Map figure

rm(list=ls())
set.seed(33)



## Simulation parametres
### S, N and SAD

s_pool = 2000
n_sim = 40000
n_mother_points = 1
max_range = 2
sampling_range = c(0.5,1.5)

### Spatial distribution (Random, medium, high)
sigma_values <- c(1, 0.1, 0.02)
sigma_names <- c("Random distribution","Medium aggregation","High aggregation")


png(filename = paste0("sampling", "S", s_pool, "_N", n_sim, "_mp", n_mother_points,".png"), width = 16, height = 5.6, unit="in", res=600)
par(mfrow = c(1,3), mex = 0.9)

lapply(1:length(sigma_values), function(sigma_i)   {
   ### Simulating the community
   comm <- mobsim::sim_thomas_community(s_pool = s_pool, n_sim = n_sim, sad_type="lnorm", fix_s_sim = TRUE, sad_coef = list("cv_abund" = 1),  sigma = sigma_values[sigma_i], mother_points = n_mother_points, xrange = c(0, max_range), yrange = c(0, max_range))
   
   
   ### Patch sizes (4 values (including very small))
   patch_areas <- c(0.0005, 0.005, 0.05, 0.5)
   patch_widths <- sqrt(patch_areas)
   
   quadrat_area <- 0.0001
   quadrat_width <- sqrt(quadrat_area)
   
   count <- 0
   maxtry <- 4000
   #### Simulating the patches and quadrats: 4 non-overlapping patches per landscape, one quadrat per patch
   ##### Patches
   
   patch_ranges <- t(sapply(patch_widths, function(patch_width) {
      patch_ranges <- c(xmin = runif(1, sampling_range[1], sampling_range[2] - patch_width), xmax = numeric(1), ymin = numeric(1), ymax = numeric(1))
      patch_ranges["xmax"] <- patch_ranges["xmin"] + patch_width
      patch_ranges["ymin"] <- runif(1, sampling_range[1], sampling_range[2] - patch_width)
      patch_ranges["ymax"] <- patch_ranges["ymin"] + patch_width
      return(patch_ranges)
   }))
   rownames(patch_ranges) <- patch_areas
   
   overlap <- any(sapply(1:length(patch_areas), function(patch_area_i){
      patch_ranges[patch_area_i,"xmin"] < patch_ranges[-patch_area_i,"xmax"] & 
         patch_ranges[patch_area_i,"xmax"] > patch_ranges[-patch_area_i,"xmin"] & 
         patch_ranges[patch_area_i,"ymin"] < patch_ranges[-patch_area_i,"ymax"] & 
         patch_ranges[patch_area_i,"ymax"] > patch_ranges[-patch_area_i,"ymin"]
   }))
   
   
   
   while(overlap && count <= maxtry){
      
      patch_ranges <- t(sapply(patch_widths, function(patch_width) {
         patch_ranges <- c(xmin = runif(1, sampling_range[1], sampling_range[2] - patch_width), xmax = numeric(1), ymin = numeric(1), ymax = numeric(1))
         patch_ranges["xmax"] <- patch_ranges["xmin"] + patch_width
         patch_ranges["ymin"] <- runif(1, sampling_range[1], sampling_range[2] - patch_width)
         patch_ranges["ymax"] <- patch_ranges["ymin"] + patch_width
         return(patch_ranges)
      }))
      rownames(patch_ranges) <- patch_areas
      
      overlap <- any(sapply(1:length(patch_areas), function(patch_area_i){
         patch_ranges[patch_area_i,"xmin"] < patch_ranges[-patch_area_i,"xmax"] & 
            patch_ranges[patch_area_i,"xmax"] > patch_ranges[-patch_area_i,"xmin"] & 
            patch_ranges[patch_area_i,"ymin"] < patch_ranges[-patch_area_i,"ymax"] & 
            patch_ranges[patch_area_i,"ymax"] > patch_ranges[-patch_area_i,"ymin"]
      }))
      
      count <- count + 1
   }
   # print(count)
   if (count > maxtry) {warning(paste0("Cannot find a sampling layout with no overlap for repetition: ", i, "
                                            Use less patches or smaller patch area."))}
   
   
   ##### Quadrat
   quadrat_xrange <- data.frame(xmin = runif(4, patch_ranges[, 1], patch_ranges[, 2] - quadrat_width))
   quadrat_xrange$xmax <- quadrat_xrange$xmin + quadrat_width
   quadrat_yrange <- data.frame(ymin = runif(4, patch_ranges[, 3], patch_ranges[, 4] - quadrat_width))
   quadrat_yrange$ymax <- quadrat_yrange$ymin + quadrat_width
   
   ## Graphical check
   plotting_range <- c(0.5, 1.5)
   colour_palette <- viridisLite::viridis(s_pool)[comm$census[,'species']] # rainbow is a basic R possibility, inferno, magma, plasma and cividis are other viridisLite possibilities (and they are all colourblind friendly palettes)
   
   plot(comm$census[,'x'], comm$census[,'y'], col = colour_palette, main = sigma_names[sigma_i], pch = 19, 
        xlim = plotting_range, ylim = plotting_range, xlab = "X", ylab = "Y", las = 1, asp = 1, xaxt = 'n', yaxt = 'n', cex = 1, xaxs = 'r', yaxs = 'r')   # plotting individual points
   # creating axes with shifted labels from 0.5-1.5 to 0-1.
   axis(side = 1,
        at = seq(from = sampling_range[1], to = sampling_range[2], by = 0.2),
        labels = seq(from = 0, to = diff(sampling_range), by = 0.2))
   axis(side = 2,
        at = seq(from = sampling_range[1], to = sampling_range[2], by = 0.2),
        labels = seq(from = 0, to = diff(sampling_range), by = 0.2), las = 1)
   graphics::rect(patch_ranges[, 1],   # 'xmin
                  patch_ranges[, 3],   # 'ymin
                  patch_ranges[, 2],   # 'xmax
                  patch_ranges[, 4],   # 'ymax
                  lwd = 2, col = grDevices::adjustcolor("white", alpha.f = 0.6))
   
   graphics::rect(quadrat_xrange[, 1],   # 'xmin
                  quadrat_yrange[, 1],   # 'ymin
                  quadrat_xrange[, 2],   # 'xmax
                  quadrat_yrange[, 2],   # 'ymax
                  lwd = 1.4, col = grDevices::adjustcolor("forestgreen", alpha.f = 0.4), border = "forestgreen")
   
})
dev.off()