forest_model <- function(t, states, parms, inputs){

  ens_members <- nrow(states)

  inputs_temp <- inputs[, 1]
  inputs_PAR <- inputs[, 2]
  inputs_doy <- inputs[, 3]

  ##Unit Conversion: umol/m2/sec to Mg/ha/timestep
  k <- 1e-6 * 12 * 1e-6 * 10000 * 86400 #mol/umol*gC/mol*Mg/g*m2/ha*sec/timestep

  ## photosynthesis
  lai <- states[, 1] * parms$SLA * 0.1  #0.1 is conversion from Mg/ha to kg/m2
  gpp <- pmax(0, k * parms$alpha * (1 - exp(-0.5 * lai)) * inputs_PAR)
  gpp[inputs_PAR < 1e-20] <- 0 ## night

  ## respiration & allocation
  ra <- gpp * parms$Ra_frac
  npp <- gpp - ra



  leaf_alloc <- npp * parms$leaf_frac
  wood_alloc <- npp * (1 - parms$leaf_frac)

  rh <- pmax(k * parms$Rbasal * states[, 3] * parms$Q10 ^ (inputs_temp / 10), 0) ## pmax ensures SOM never goes negative

  ## turnover


  litterfall <- states[ , 1] * (parms$litterfall_rate * (365/ (params$litterfall_length)))
  litterfall[!(inputs_doy > params$litterfall_start & inputs_doy[1] < (params$litterfall_start + params$litterfall_length))] <- 0.0

  mortality <- states[ , 2] * parms$mortality

  #Change in states
  dleaves_dt <- leaf_alloc  - litterfall
  dwood_dt <- wood_alloc  - mortality
  dSOM_dt <- litterfall + mortality - rh

  #
  states[, 1] <- states[, 1] + dleaves_dt
  states[, 2] <- states[, 2] + dwood_dt
  states[, 3] <- states[, 3] + dSOM_dt

  ## update states
  states[, 1] <- pmax(rnorm(ens_members, states[, 1] , parms$sigma.leaf), 0)
  states[, 2] <- pmax(rnorm(ens_members, states[, 2], parms$sigma.wood), 0)
  states[, 3] <- pmax(rnorm(ens_members, states[, 3], parms$sigma.soil), 0)

  #Extra variables
  lai <- states[, 1]  * parms$SLA * 0.1
  nee <- ra + rh - gpp

  return(cbind(state1 = states[, 1],
               state2 = states[, 2],
               state3 = states[, 3],
               lai = lai,
               gpp = gpp ,
               nee = nee,
               ra =  ra,
               npp_w = wood_alloc,
               npp_l = leaf_alloc,
               rh = rh,
               litterfall = litterfall,
               mortality = mortality))

}
