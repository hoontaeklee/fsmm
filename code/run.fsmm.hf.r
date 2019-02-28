## ----------------------------------------------------------------------- ---- 
## initial setting ------------------------------------------------------- ----
rm(list = ls())
# install.packages("suncalc")

setwd("D:/fmc_model/fsmm-master/")

## ----------------------------------------------------------------------- ---- 
## constants ----

## albedo (Table D-2 Physical Hydrology 2nd ed. Lawrence Dingman)
albedo_hfo = 0.2 # Grass - typical field
albedo_hfi = 0.1 # Red pine forest

## emissivity (Table D-1 Physical Hydrology 2nd ed. Lawrence Dingman)
emissivity_hfo = 0.95 # Grass - typical field
emissivity_hfi = 0.97 # coniferous forest

## conversion factor (MJ/m2 --> W/m2)
fac_mj_w = 1 / 0.0036

## data load ----
df_hfo = read.csv("./df_test/df_hfo.csv",
                  header = TRUE, stringsAsFactors = FALSE)
df_hfo$date = as.POSIXct(df_hfo$date,
                         format = "%Y-%m-%d %H:%M:%S",
                         tz = "UTC")

df_hfo$press = as.numeric(gsub(",", "", df_hfo$press))

df_hfi = read.csv("./df_test/df_hfi.csv",
                  header = TRUE, stringsAsFactors = FALSE)
df_hfi$date = as.POSIXct(df_hfi$date,
                         format = "%Y-%m-%d %H:%M:%S",
                         tz = "UTC")
df_hfi$press = as.numeric(gsub(",", "", df_hfi$press))

## pressure at the sea level
p.sea = read.csv("./df_test/ASOS_108_Psea_HR_2018080100_2019022523.csv",
                 header = TRUE, stringsAsFactors = FALSE)
p.sea$date = as.POSIXct(p.sea$date,
                        format = "%Y-%m-%d %H:%M")
seq_date_p.sea = 
  data.frame("date" = seq.POSIXt(from = p.sea$date[1],
                                 to = p.sea$date[dim(p.sea)[1]],
                                 by = "1 hour"))

## expandinig the NAs
p.sea =
  merge(seq_date_p.sea, p.sea, 
        by = "date", all.x = TRUE)

## omit dates not in observation
p.sea = p.sea[p.sea$date %in% df_hfo$date, ]


## total shortwave radiation 
total_swdn = read.csv("./df_test/ASOS_108_Rad_HR_2018080100_2019022523.csv",
                      header = TRUE, stringsAsFactors = FALSE)
total_swdn$k.down = total_swdn$k.down * fac_mj_w
total_swdn$date = as.POSIXct(total_swdn$date,
                             format = "%Y-%m-%d %H:%M")
seq_date_swdn = 
  data.frame("date" = seq.POSIXt(from = total_swdn$date[1],
                                 to = total_swdn$date[dim(total_swdn)[1]],
                                 by = "1 hour"))

## expandinig the NAs
total_swdn =
  merge(seq_date_swdn, total_swdn, 
        by = "date", all.x = TRUE)

## swdn NA --> 0 (no measurement at night)
total_swdn[is.na(total_swdn$k.down) == TRUE, "k.down"] = 0

## omit dates not in observation
total_swdn = total_swdn[total_swdn$date %in% df_hfo$date, ]


## ----------------------------------------------------------------------- ---- 
## make input data (direct swdn, diffuse swdn, lwdn)---------------------- ----
## ----------------------------------------------------------------------- ---- 
## solar altitude and azimuth angle ----
df_hfo$Alt_s = 
  oce::sunAngle(t = df_hfo$date, 
                longitude = 127.042806,
                latitude = 37.593806)$altitude * pi / 180
df_hfo$Azi_s = 
  oce::sunAngle(t = df_hfo$date,
                longitude = 127.042806,
                latitude = 37.593806)$azimuth * pi / 180

df_hfi$Alt_s = 
  oce::sunAngle(t = df_hfi$date,
                longitude = 127.045197, 
                latitude = 37.595600)$altitude * pi / 180
df_hfi$Azi_s = 
  oce::sunAngle(t = df_hfi$date,
                longitude = 127.045197,
                latitude = 37.595600)$azimuth * pi / 180

## downward shortwave radiation (direct + diffuse) ----

## a function computing the maximum downward shortwave radiation
calc.k.down.max = function(p, # pressure at the site
                           p.sea, # pressure at the sea-level
                           phi, # solar elevation angle
                           tau = 0.75, # atmospheric transmissivity
                           k.solar = 1367){ # solar constant (Wm-2)
  k.down.max = k.solar * tau ^ (p / p.sea) * sin(phi)
  
  k.down.max[k.down.max < 0] = 0
  
  return(k.down.max)
} 

## a function computing the fraction of diffuse downward shortwave radiation
calc.frac.k.down.diff = function(n){ 
  # n: k.down / k.down.max
  # frac.k.down.diff: k.down.diff / k.down
  
  frac.k.down.diff = n
  
  idx_under_0.22 = n <= 0.22
  frac.k.down.diff[idx_under_0.22] =
    1.0 - 0.09 * frac.k.down.diff[idx_under_0.22]
  
  idx_over_0.22 = n > 0.22 & n <= 0.22
  frac.k.down.diff[idx_over_0.22] =
    0.951 - 0.1604 * frac.k.down.diff[idx_over_0.22] + 
    4.388 * frac.k.down.diff[idx_over_0.22] ^ 2 -
    16.638 * frac.k.down.diff[idx_over_0.22] ^ 3 +
    12.336 * frac.k.down.diff[idx_over_0.22] ^ 4
  
  idx_over_0.80 = n > 0.80
  frac.k.down.diff[idx_over_0.80] = 0.165

  frac.k.down.diff[frac.k.down.diff > 1] = 1
  frac.k.down.diff[frac.k.down.diff < 0] = 0
  
  return(frac.k.down.diff)
}


## main calculation
df_hfo$k.down.max = calc.k.down.max(p = df_hfo$press,
                                    p.sea = p.sea$p.sea,
                                    phi = df_hfo$Alt_s)
df_hfo$n = total_swdn$k.down / df_hfo$k.down.max
df_hfo$n[is.na(df_hfo$n) | is.infinite(df_hfo$n)] = 0

df_hfo$frac.k.down.diff = calc.frac.k.down.diff(df_hfo$n)

df_hfo$k.down.dir = total_swdn$k.down * (1 - df_hfo$frac.k.down.diff)
df_hfo$k.down.diff = total_swdn$k.down * df_hfo$frac.k.down.diff

# df_hfo$k.sw = 
#   conv_factor * total_swdn$k.down / albedo_hfo  
# 
# df_hfi$k.sw = 
#   conv_factor * total_swdn$k.down / albedo_hfi  

## downward longwave radiation ----
calc.ea = function(ta, rh){
  # reference: plants and microclimate 2nd ed.
  # es: saturated vapor pressure (hPa)
  # ea: vapor pressure (hPa)
  # ta: Celsius air temperature
  # rh: relative humidity (-)
  
  a = 613.75
  b = 17.502
  c = 240.97
  
  es = a * exp(b * ta / (ta + c)) * 0.01
  ea = es * rh
  
  return(ea)
}

calc.emis.sky = function(ea, ta){
  # w: the precipitable water content (cmm)
  # M: molecular mass of water (0.0180 kg mol-1) 
  # R: gas constant (8.314 * 10e-3 kPa m3 mol-1 K-1)
  # ea: vapor pressure (hPa)
  # ta: Celsius air temperature
  
  M = 0.0180
  R = 8.314 * 10e-3
  w = 46.5 * ea / (ta + 273.15)
  
  emis.sky = 1 - (1 + w) * exp(-(1.2 + 3.0 * w) ^ 0.5)
  
  return(emis.sky)
}


calc.emis.atm = function(emis.sky, n){
  # emis.sky: emissivity of clear sky
  # n: clearness index (k.down / k.down.max)
  # beta: the cloudiness factor
  
  beta = 0.26
  
  m = n
  
  m[n > 0.2] = 1 - n
  m[n <= 0.2] = 1
  
  emis.atm = (1 + beta * m) * emis.sky
  emis.atm[emis.atm > 1] = 1
  return(emis.atm)
}

calc.l.dn = function(s, 
                     emis.atm,
                     emis.canopy,
                     ta){
  # l.dn: downward longwave radiation (W m-2)
  # s: sky view factor
  # emis.atm = atmospheric emissivity 
  # emis.canopy = canopy emissivity
  # SB.con = stephan-Boltzmann constant (5.67 * 1e-8 W m-2 K-4)
  # ta: Celsius air temperature
  
  SB.con = 5.67 * 1e-8
  l.dn = (s * emis.atm + (1 - s) * emis.canopy) * SB.con * (ta + 273.15) ^ 4
  
  return(l.dn)
}

df_hfo$ea = calc.ea(df_hfo$ta_2_avg,
                    df_hfo$rh_2_avg)
df_hfo$emis.sky = calc.emis.sky(df_hfo$ea,
                                df_hfo$ta_2_avg)
df_hfo$emis.atm = calc.emis.atm(df_hfo$emis.sky,
                                df_hfo$n)
df_hfo$l.dn = calc.l.dn(s = 1.0,
                        emis.atm = df_hfo$emis.atm,
                        emis.canopy = emissivity_hfo,
                        ta = df_hfo$ta_2_avg)

## ----------------------------------------------------------------------- ----
## run the fsmm model  ----
input_fsmm_hfo =
  data.frame("Alt_s" = df_hfo$Alt_s,
             "Azi_s" = df_hfo$Azi_s,
             "k.down.dir" = df_hfo$k.down.dir,
             "k.down.diff" = df_hfo$k.down.diff,
             "temp" = df_hfo$ta_2_avg,
             "L_dn" = df_hfo$l.dn,
             "rh" = df_hfo$rh_2_avg,
             "wind.speed" = df_hfo$ws_2_avg,
             "precip" = df_hfo$rain)
fsmm.model.input = 
  read.csv(file.path(src.dir,
                     "data/fsmm.example.forcing.vars.csv"))

summary(fsmm.model.input$k.down.diff)
summary(input_fsmm_hfo$k.down.diff)

summary(fsmm.model.input$k.down.dir)
summary(input_fsmm_hfo$k.down.dir)


summary(total_swdn$k.down)
c(6.7, -10, 0.1e-10, 0.35, 0.22)
c(4.4,-10,3e-10,0.3,0.1)
c(5.164862,-19.84168,4.176724e-10,0.3,0.16542898) # NSE.log = -0.3816146

output <- 
  run.fsmm(fsmm.model.input = input_fsmm_hfo,
             fsmm.model.pars = c(4.505667e+00,-1.000000e+01,4.999938e-10,0.3,0.9),
           ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
           RHO_S = 400,
           # kg m-3 : Stick density  (taken from Nelson, 2000))
           r = 0.0125,
           # radius of stick (m)
           L = 0.13,
           # length of stick (m)
           m_i = 0.1,
           CALC_SHAPE_FACTOR = T)
df_hfo$fmc_fsmm = output$fuel.mois.mod
plot(df_hfo$fmc_fs3, df_hfo$fmc_fsmm * 100,
     xlim = c(0, 35), ylim = c(0, 35))
abline(0, 1)

points(df_hfo$fmc_fs3, type = "l")
plot(df_hfo$fmc_fsmm * 100, col = "red", type = "l")
summary(lm(df_hfo$fmc_fsmm * 100 ~ df_hfo$fmc_fs3))
nse = 1 - (sum((log(df_hfo$fmc_fsmm * 100) - log(df_hfo$fmc_fs3)) ^ 2) / sum((log(df_hfo$fmc_fsmm * 100) - mean(log(df_hfo$fmc_fs3))) ^ 2))
## ----------------------------------------------------------------------- ----
## optimize with a particle-swarm optimization routine ----

calc.nse = function(input.fsmm = input.fsmm,
                    params = params,
                    ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                    RHO_S = RHO_S,
                    # kg m-3 : Stick density  (taken from Nelson, 2000))
                    r = r,
                    # radius of stick (m)
                    L = L,
                    # length of stick (m)
                    m_i = m_i,
                    CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR,
                    obs){
  -1 + 
    (sum((run.fsmm(fsmm.model.input = input.fsmm,
                   fsmm.model.pars = params,
                   ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                   RHO_S = RHO_S,
                   # kg m-3 : Stick density  (taken from Nelson, 2000))
                   r = r,
                   # radius of stick (m)
                   L = L,
                   # length of stick (m)
                   m_i = m_i,
                   CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR)$fuel.mois.mod * 100 - obs) ^ 2) /
       sum((run.fsmm(fsmm.model.input = input.fsmm,
                     fsmm.model.pars = params,
                     ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                     RHO_S = RHO_S,
                     # kg m-3 : Stick density  (taken from Nelson, 2000))
                     r = r,
                     # radius of stick (m)
                     L = L,
                     # length of stick (m)
                     m_i = m_i,
                     CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR)$fuel.mois.mod * 100 - mean(obs)) ^ 2))
}

calc.nse.log = function(input.fsmm = input.fsmm,
                        params = params,
                        ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                        RHO_S = RHO_S,
                        # kg m-3 : Stick density  (taken from Nelson, 2000))
                        r = r,
                        # radius of stick (m)
                        L = L,
                        # length of stick (m)
                        m_i = m_i,
                        CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR,
                        obs){
  -1 + 
    (sum((log(run.fsmm(fsmm.model.input = input.fsmm,
                       fsmm.model.pars = params,
                       ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                       RHO_S = RHO_S,
                       # kg m-3 : Stick density  (taken from Nelson, 2000))
                       r = r,
                       # radius of stick (m)
                       L = L,
                       # length of stick (m)
                       m_i = m_i,
                       CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR)$fuel.mois.mod * 100) - obs) ^ 2) /
       sum((log(run.fsmm(fsmm.model.input = input.fsmm,
                         fsmm.model.pars = params,
                         ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                         RHO_S = RHO_S,
                         # kg m-3 : Stick density  (taken from Nelson, 2000))
                         r = r,
                         # radius of stick (m)
                         L = L,
                         # length of stick (m)
                         m_i = m_i,
                         CALC_SHAPE_FACTOR = CALC_SHAPE_FACTOR)$fuel.mois.mod * 100) - mean(obs)) ^ 2))
}

## optim against nse
pso::psoptim(par = rep(NA, 5), 
             fn = calc.nse,
             input.fsmm = input_fsmm_hfo,
             RHO_S = 400,
             r = 0.0125,
             L = 0.13,
             m_i = 0.1,
             CALC_SHAPE_FACTOR = T,
             obs = df_hfo$fmc_fs3,
             lower = c(4.4, -22, 0.1e-10, 0.3, 0.05),
             upper = c(6.7, -10, 5.0e-10, 0.3, 0.9))

## optim against nse.log
nse.log = -10
cnt_consecutive = 0
cnt_iter = 0
mat = matrix(NA, ncol = 7, nrow = 2)
colnames(mat) = c("index", "nse.log", "par1", "par2", "par3", "par4", "par5")
mat = as.data.frame(mat)
write.csv(mat, "D:/psoptim_hfo.csv", row.names = FALSE)
while(cnt_consecutive < 20){
  cnt_iter = cnt_iter + 1
  
  run = 
    pso::psoptim(par = rep(NA, 5), 
                 fn = calc.nse.log,
                 input.fsmm = input_fsmm_hfo,
                 RHO_S = 400,
                 r = 0.0125,
                 L = 0.13,
                 m_i = 0.1,
                 CALC_SHAPE_FACTOR = T,
                 obs = log(df_hfo$fmc_fs3),
                 lower = c(4.4, -22, 0.1e-10, 0.3, 0.05),
                 upper = c(6.7, -10, 5.0e-10, 0.3, 0.9),
                 control = list(abstol = -1))
  
  improvement = run$value - nse.log
  
  mat[cnt_iter, ] = c(cnt_iter, run$value, run$par)
  
    if(improvement >= 1e-4){
    nse.log = nse.log + improvement
    cnt_consecutive = 0
  } else{
    cnt_consecutive = cnt_consecutive + 1
  }

  
  print(paste0("current nse = ", nse.log))
  print(paste0("cnt_iter = ", cnt_iter))
  print(paste0("cnt_consecutive = ", cnt_consecutive))

  }

## ----------------------------------------------------------------------- ----

