# Review - Simulation


## Roadmap

### Simulation parametres
#### S, N and SAD

#### Spatial distribution (Random, medium, high aggregation)
##### Per species ranges (not implemnented yet)

##### Aggregation parametres

#### Patch sizes (4 values (including very small))

### Running simulations
#### Simulating the patch

#### sampling inside the patch

#### Saving N, S and S_PIE for each simulation

### Products
#### Maps (including some with sample sizes)

#### Distribution of lm coefficients

#### Saved simulations summary values

####################################################################

rm(list=ls())
set.seed(33)

## Simulation parametres
### S, N and SAD

s_pool = 2000
n_sim = 40000
n_mother_points = 1
max_range = 2
sampling_range = c(0.5,1.5)

metrics <- c("N", "S", "S_PIE")
nrep = 2000


### Spatial distribution (Random, medium, high)
#### Aggregation parametres
sigma_values <- c(1, 0.1, 0.02)
#### Computing and storing community tables
communities <- array(NA, dim = c(n_sim, 3, length(sigma_values), nrep), 
                     dimnames = list(c(), c('x','y','species'), paste0("sigma_", sigma_values), c())
                     )
for (sigma_i in 1:length(sigma_values))   {
   for(i in 1:nrep)   {
      communities[,, sigma_i, i] <- unlist(
			mobsim::sim_thomas_community(s_pool = s_pool, n_sim = n_sim, sad_type="lnorm", 
			                             fix_s_sim = TRUE, sad_coef = list("cv_abund" = 1),  
			                             sigma = sigma_values[sigma_i], mother_points = n_mother_points, 
			                             xrange = c(0, max_range), yrange = c(0, max_range))$census
		)
   }
}


##### map of the community

# svg(filename = paste0(path2wd, "intermediate_results/", "7_mapS", s_pool, "_N", n_sim, "_mp", n_mother_points,".svg"), width = 20, height = 8)
# par(mfrow=c(1,3))
# lapply(1:length(sigma_values), function(sigma_i) {
# plot(communities[,'x',sigma_i, 1], communities[,'y',sigma_i, 1], 
#      main = paste("sigma =", sigma_values[sigma_i]), 
#      cex = 0.8, las = 1, asp = 1, col = communities[,'species',sigma_i, 1], pch = 19)
# })
# dev.off()

### Patch sizes (4 values (including very small))
patch_areas <- c(0.0005, 0.005, 0.05, 0.5)
patch_widths <- sqrt(patch_areas)

quadrat_area <- 0.0001
quadrat_width <- sqrt(quadrat_area)


## Running simulations

res <- array(data = NA,
             dim = list(nrep, length(patch_areas), length(sigma_values), length(metrics)),
             dimnames = list(1:nrep, paste0("patch_area_", patch_areas), paste0("sigma_", sigma_values), metrics))


for(i in 1:nrep)   {
   for (sigma_i in 1:length(sigma_values))   {
 
      comm <- communities[,, sigma_i, i]
      
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
      if(FALSE) {
         plot(comm[,'y'] ~ comm[,'x'], col = comm[,'species'], main = "Community distribution", las = 1, asp = T, pch=20)
         graphics::rect(patch_ranges[, 1],
                        patch_ranges[, 3],
                        patch_ranges[, 2],
                        patch_ranges[, 4],
                        lwd = 2, col = grDevices::adjustcolor("white", alpha.f = 0.6))
         
         graphics::rect(quadrat_xrange[, 1],
                        quadrat_yrange[, 1],
                        quadrat_xrange[, 2],
                        quadrat_yrange[, 2],
                        lwd = 1.4, col = grDevices::adjustcolor("forestgreen", alpha.f = 0.4), border = "forestgreen")
      }
      
      #### sampling inside the quadrat
      quadrat_comms <- lapply(1:length(patch_widths), function(patch_width_i){
         comm[comm[,'x'] >= quadrat_xrange[patch_width_i, 1] & comm[,'x'] <= quadrat_xrange[patch_width_i, 2] & comm[,'y'] >= quadrat_yrange[patch_width_i, 1] & comm[,'y'] >= quadrat_yrange[patch_width_i, 2],]
      })  
      
      ### Saving N, S and S_PIE for each simulation
      res[i, , sigma_i, ] <-  c(t(sapply(quadrat_comms, function(quadrat_comm){
         c(N = nrow(quadrat_comm),
           S = length(unique(quadrat_comm[,'species'])),
           S_PIE = mobr::calc_PIE(ENS = TRUE, table(quadrat_comm[,'species'])))
      })))
   }
}

# It seems this is not necessary right?
# save(res, file = paste0(path2wd, "intermediate_results/7_raw_results", "S", s_pool, "_N", n_sim, "_mp", n_mother_points, "_nrep", nrep))

# saving results in long format
restable <- as.data.frame.table(res)
colnames(restable) <- c("rep","patch_area","sigma","metric","value")

restable$patch_area <- as.numeric(gsub(restable$patch_area, pattern = "patch_area_", replacement = ""))
restable$sigma <- as.numeric(gsub(restable$sigma, pattern = "sigma_", replacement = ""))

write.table(restable, 
            file = paste0(path2wd, "intermediate_results/7_results", "S", s_pool, "_N", n_sim, "_mp", n_mother_points, "_nrep", nrep, ".csv"), 
            row.names = F, sep = ',', dec='.')

# I think the only bug here was a missing slash ;o)
write.table(restable, file = paste0(path2wd, "/intermediate_results/7_results", "S", s_pool, "_N", n_sim, "_mp", n_mother_points, "_nrep", nrep, ".csv"), row.names = F, sep = ',', dec='.')


