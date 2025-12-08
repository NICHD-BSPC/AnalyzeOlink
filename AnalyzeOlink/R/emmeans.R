# lmrob <-> emmeans compatiblity wrappers
recover_data.lmrob <- function(object, ...)
  emmeans:::recover_data.lm(object, ...)

emm_basis.lmrob <- function(object, trms, xlev, grid, ...)
  emmeans:::emm_basis.lm(object, trms, xlev, grid, ...)
