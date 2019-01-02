# code to examine pair-wise comparisons for the two-way interactions:
# Questions: 
# 1) do the responses to fragment size differ between taxa? 
# 2) does the response to fragment size depend on the harshness of the matrix?
# 3) does the response to fragment size depend on the time since fragmentation?
# 4) does the response to fragment size depend on the fragment vegetation?


# load the model fits and the brms library

load('~/Dropbox/Habitat loss meta-analysis/analysis/two_interactions_brms_fits.Rdata')

library(brms)


# ---TAXA: are the relationships with fragment size the same for all taxa? Specify all pairwise comparisons-----
taxa_hyp <- c('c.lfs = c.lfs:taxabirds',
       'c.lfs = c.lfs:taxainvertebrates',
       'c.lfs = c.lfs:taxamammals',
       'c.lfs = c.lfs:taxaplants',
       'c.lfs:taxabirds = c.lfs:taxainvertebrates',
       'c.lfs:taxabirds = c.lfs:taxamammals',
       'c.lfs:taxabirds = c.lfs:taxaplants',
       'c.lfs:taxainvertebrates = c.lfs:taxamammals',
       'c.lfs:taxainvertebrates = c.lfs:taxaplants',
       'c.lfs:taxamammals = c.lfs:taxaplants')

N_fS_taxa_hypoth <- hypothesis(lNstd_fS_taxa,
                               hypothesis = taxa_hyp)

S_std_fS_taxa_hypoth <- hypothesis(lS_std_fS_taxa,
                               hypothesis = taxa_hyp)

S_PIE_fS_taxa_hypoth <- hypothesis(lS_PIE_fS_taxa,
                               hypothesis = taxa_hyp)

S_n_fS_taxa_hypoth <- hypothesis(lSn_fS_taxa,
                                   hypothesis = taxa_hyp)


# ---matrix harshness--------
matrix_hyp <- c('c.lfs = c.lfs:matrix.categorymediumfilter',
                'c.lfs = c.lfs:matrix.categoryharshfilter')

N_fS_matrix_hypoth <- hypothesis(lNstd_fS_matrix,
                               hypothesis = matrix_hyp)

S_std_fS_matrix_hypoth <- hypothesis(lS_std_fS_matrix,
                                   hypothesis = matrix_hyp)

S_PIE_fS_matrix_hypoth <- hypothesis(lS_PIE_fS_matrix,
                                   hypothesis = matrix_hyp)

S_n_fS_matrix_hypoth <- hypothesis(lSn_fS_matrix,
                                 hypothesis = matrix_hyp)



# ---time since fragmentation--------
time_hyp <- c('c.lfs = c.lfs:time.since.fragmentationintermediate20M100years',
                'c.lfs = c.lfs:time.since.fragmentationlong100Pyears')

N_fS_time_hypoth <- hypothesis(lNstd_fS_time,
                                 hypothesis = time_hyp)

S_std_fS_time_hypoth <- hypothesis(lS_std_fS_time,
                                     hypothesis = time_hyp)

S_PIE_fS_time_hypoth <- hypothesis(lS_PIE_fS_time,
                                     hypothesis = time_hyp)

S_n_fS_time_hypoth <- hypothesis(lSn_fS_time,
                                   hypothesis = time_hyp)




# ---fragment vegetation--------
veg_hyp <- c('c.lfs = c.lfs:veg.fragmentgrassland',
              'c.lfs = c.lfs:veg.fragmentshrublandDsteppe')

N_fS_veg_hypoth <- hypothesis(lNstd_fS_veg,
                               hypothesis = veg_hyp)

S_std_fS_veg_hypoth <- hypothesis(lS_std_fS_veg,
                                   hypothesis = veg_hyp)

S_PIE_fS_veg_hypoth <- hypothesis(lS_PIE_fS_veg,
                                   hypothesis = veg_hyp)

S_n_fS_veg_hypoth <- hypothesis(lSn_fS_veg,
                                 hypothesis = veg_hyp)

