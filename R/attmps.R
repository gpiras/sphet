#' #' #' Generate impacts for spreg lag and sarar models
#' #' #' 
#' #' #' 
#' #' #' @param object An object of class sphet
#' #' #' @param ... Additional arguments to be passed 
#' #' #' 
#' #' #' 
#' #' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' #' 
#' #' #' @examples 
#' #' #' data(columbus, package="spdep")
#' #' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' #' summary(res)
#' #' #' effects <- impacts(res)
#' #' #' summary(effects)
#' #' impacts <- function(object, ...){
#' #'   UseMethod("impacts", object)
#' #' }
#' #' #' describedIn impacts
#' #' #' @param object An object of class sphet
#' #' #' @param ... Additional arguments to be passed 
#' #' #' 
#' #' #' 
#' #' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' #' 
#' #' #' @examples 
#' #' #' data(columbus, package="spdep")
#' #' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' #' summary(res)
#' #' #' effects <- impacts(res)
#' #' #' summary(effects)
#' #' 
#' #' 
#' 
#' #    rm(list = ls())
#' #    ############################################
#' #    ##### spreg "lag"
#' #    ############################################
#' #
#' #    source("/home/gpiras/R dev lib/sphet/R/impcts.R")
#' #    #  source("/home/gpiras/R dev lib/spatialreg/R/impacts.R")
#' #
#'      data(columbus, package =  "spdep")
#' 
#'      listw <- spdep::nb2listw(col.gal.nb)
#'      wi <- spatialreg::eigenw(listw)
#' # #### lag model
#'     res <- sphet::spreg(CRIME~INC +HOVAL, data=columbus ,
#'                         listw= listw, model = "lag" , Durbin=F)
#'    #listw
#'   # set.seed("1234")
#'   # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#'   set.seed("1234")
#'    summary(impacts(res, listw = listw, R = 10))
#'    summary(eff)
#' #
#' #    #eigenvalues
#' #    set.seed("1234")
#' #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T, weg = wi))
#'     set.seed("1234")
#'     eff <- spatialreg::impacts(res,  R = 10, evalues = wi)
#'     summary(eff)
#' #
#' #
#' #    #trace
#'     wm <- spdep::listw2mat(listw)
#'     ws <- as(wm, "sparseMatrix")
#'     traccia <- spatialreg::trW(ws)
#' #    set.seed("123444")
#' #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = F, trW = traccia))
#'     set.seed("123444")
#'     eff <- spatialreg::impacts(res,  R = 10, tr = traccia)
#'     summary(eff)
#' #
# ######spatial Durbin Model
# ##### Durbin = TRUE
#    #    rm(list=ls())
  # data(boston, package = "spData")
  # listw <- spdep::nb2listw(boston.soi)
  # res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "lag", Durbin = TRUE)
#  set.seed("45")
#  summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#  set.seed("45")
  # eff <- spatialreg::impacts(res, listw = listw, R = 10)
#  summary(eff)
# #
#  ##### Durbin = ~RM (only one of the x)
#
  # res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "lag", Durbin = ~RM)
#   set.seed("45")
#    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    set.seed("45")
   #  eff <- spatialreg::impacts(res, listw = listw, R = 10)
   # summary(eff)
#
#
#    ##### Durbin = ~CRIM +RM (only one of the x)
#
    # res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "lag", Durbin = ~CRIM + RM)
#    set.seed("45")
#    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    set.seed("45")
#     eff <- spatialreg::impacts(res, listw = listw, R = 10)
#     summary(eff)
# #
#
#    ##### Durbin = ~CRIM +RM (only two of the x)
#
#     res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "lag", Durbin = ~RM + CRIM )
# #    set.seed("45")
# #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
# #    set.seed("45")
#     eff <- spatialreg::impacts(res, listw = listw, R = 10)
#     summary(eff)
# #
#
#
#    ##### Durbin = ~CRIM +RM (a wx which is not in x)
#
    # res <- spreg(CMEDV ~AGE+ RM, data = boston.c, listw = listw,  model = "lag", Durbin = ~RM + CRIM )
#    set.seed("45")
#    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    set.seed("45")
    # eff <- spatialreg::impacts(res, listw = listw, R = 10)
    # summary(eff)
#
#
#
#
   ############################################
   ##### spreg "sarar"
#    ############################################
# 
# # rm(list=ls())
# #    source("/home/gpiras/R dev lib/sphet/R/impcts.R")
#     data(columbus, package =  "spdep")
#     listw <- spdep::nb2listw(col.gal.nb)
#     wi <- spatialreg::eigenw(listw)
# #    #### lag model
#     res <- sphet::spreg(CRIME~INC +HOVAL, data=columbus ,
#                         listw= listw, model = "sarar" , Durbin=F)
# #    #listw
# #    set.seed("1234")
# #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
# #    set.seed("1234")
#     eff <- spatialreg::impacts(res, listw = listw, R = 10)
#     summary(eff)
# #
# #    #eigenvalues
# #    set.seed("1234")
# #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T, weg = wi))
# #    set.seed("1234")
#     eff <- spatialreg::impacts(res,  R = 10, evalues = wi)
#     summary(eff)
# #
# #
# #    #trace
#     wm <- spdep::listw2mat(listw)
#     ws <- as(wm, "sparseMatrix")
#     traccia <- spatialreg::trW(ws)
# #    set.seed("123444")
# #    summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = F, trW = traccia))
# #    set.seed("123444")
#     eff <- spatialreg::impacts(res,  R = 10, tr = traccia)
#     summary(eff)
# #
#     ######spatial Durbin Model
#    ##### Durbin = TRUE
#    #    rm(list=ls())
#  data(boston, package = "spData")
#     listw <- spdep::nb2listw(boston.soi)
#     res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "sarar", Durbin = TRUE)
#    # set.seed("45")
#    # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    # set.seed("45")
#     eff <- spatialreg::impacts(res, listw = listw, R = 10)
#     summary(eff)
#    # #
#    # ##### Durbin = ~RM (only one of the x)
#    #
#     res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "sarar", Durbin = ~RM)
#    
#     # set.seed("45")
#    # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    # set.seed("45")
#     eff <- sphet::impacts.gstsls(res, listw = listw, R = 10)
#     summary(eff)
#    #
#    # #
#    # ##### Durbin = ~CRIM +RM (only one of the x)
#    #
#     res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "sarar", Durbin = ~CRIM + RM)
#    # set.seed("45")
#    # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
#    # set.seed("45")
#     eff <- sphet::impacts(res, listw = listw, R = 10)
#     summary(eff)
#    #
#    #   # ##### Durbin = ~CRIM +RM (only two of the x)   
#   # res <- spreg(CMEDV ~AGE+ CRIM +RM, data = boston.c, listw = listw,  model = "sarar", Durbin = ~RM + CRIM )   
    # set.seed("45")
    # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))
    # set.seed("45")   
    # eff <- impacts.gstsls(res, listw = listw, R = 10)   
    # summary(eff)
    #   #   #   # ##### Durbin = ~CRIM +RM (a wx which is not in x)
    #   # res <- spreg(CMEDV ~AGE+ RM, data = boston.c, listw = listw,  model = "sarar", Durbin = ~RM + CRIM )  
    # set.seed("45")   
    # summary(sphet::impacts(res, n_mvn = 10, inference = T, exact = T))   
    # set.seed("45")
    # eff <- impacts.gstsls(res, listw = listw, R = 10)
    # summary(eff)