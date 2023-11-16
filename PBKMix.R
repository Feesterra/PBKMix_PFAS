### Mixture PBK model for PFAS for zebrafish female ###
# 11/2023
# E. Golosovskaia, S. Örn, L. Ahrens, I. Chelcea, P. L Andersson


##############################
### Variable -(Notation)- Units
# Amounts:        -(A_x) - microg 
# Volumes:      -(V_x) - mL
# Time:            -(t)   - d
# Flows:           -(F_x) - mL/d
# Concentrations: -(C_x) - microg/mL
# Masses:         -(BW)  - g
# Temperature: -(TC_c)  - Celsius
# Temperature: -(TC_k)  - Kelvin
# Density of each tissue is considered equal to 1
###############################


#### Arrhenius temperatures function ### as presented by Grech et al. 2019
KT <- function(T, TR, TA){ exp ( (TA / TR) - (TA / T) ) }

#===============================================
# Model
#===============================================

ZF.model.c <- function(t,initial_v, phys_exp_parms, chem_parms) {
  with(as.list(phys_exp_parms),{
    
    # Declaring state vars
    y<-unname(initial_v)
    
    V_egg = y[1] # N1
    V_urine = y[2] # N2
    
    k=0
    A_admin_gill = y[(3+n*k):((3+n-1)+n*k)]  # N3:4
    k=k+1
    A_lumen = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_excr_gill = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_urine = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_urine_cum = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_feces = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_excr_water = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_bile = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_blood = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_liv = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_egg = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_egg_cum = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_gon = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_fat = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_git = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_brain = y[(3+n*k):((3+n-1)+n*k)]  
    k=k+1
    A_kidney = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_skin = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_rp = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_pp = y[(3+n*k):((3+n-1)+n*k)] 
    k=k+1
    A_water = y[(3+n*k):((3+n-1)+n*k)] 
    
    # Defining chemical parameters
    chem_p = unname(chem_parms)
    k=0
    i=chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    pKa = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    logKow = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    MW = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Dmw = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Dpw = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    pKd = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pegggon = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Plivb = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Pgonb =  chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Pbb = chem_p[(1+n*k):((1+n-1)+n*k)] 
    k=k+1
    Prpb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pppb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pgitb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pfatb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pskinb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    Pkidb = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    
    #Oral absorption
    Ku       = chem_p[(1+n*k):((1+n-1)+n*k)]   # Diffusion coefficient
    k=k+1
    frac_abs = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    
    # Excretion
    Cl_bile  = chem_p[(1+n*k):((1+n-1)+n*k)] # ml/d
    k=k+1
    Ke_urine = chem_p[(1+n*k):((1+n-1)+n*k)] # 1/d
    k=k+1
    K_BG     = chem_p[(1+n*k):((1+n-1)+n*k)] # 1/d
    k=k+1
    Ke_feces = chem_p[(1+n*k):((1+n-1)+n*k)]  # 1/d estiamted from Nichols  et al. 2004 0.83
    
    k=k+1
    PL = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    blood_bound = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    WaterExposure  = chem_p[(1+n*k):((1+n-1)+n*k)]
    k=k+1
    
    ## Scalar parameter calculation
    
    dV_egg = Egggrowth
    
    PS       = pi^(1/3)*(6*V_one_egg)^(2/3)*(V_egg/V_one_egg)

    #Volumes (ml)
    V_blood = sc_blood * BW
    V_liv    = sc_liv    * BW
    V_gon    = sc_gon    * BW
    V_fat    = sc_fat    * BW
    V_git    = sc_git    * BW
    V_brain  = sc_brain  * BW
    V_kidney = sc_kidney * BW
    V_skin   = sc_skin   * BW
    V_rp     = sc_rp     * BW
    V_pp     = BW-(V_blood+V_liv+V_gon+V_fat+V_git+V_brain+V_kidney+V_skin+V_rp)
    
    V_total = V_blood+V_liv+V_gon+V_fat+V_git+V_brain+V_kidney+V_skin+V_rp+V_pp
    V_bal = BW- V_total
    
    #Parameter equations
    TC_k = TC_c + 273.15 # (degree K)
    Fcard = (F_card_ref * KT(T=TC_k , TR=TR_Fcard,TA=TA) * (BW/Bw_Fcard_ref)^(-0.1))*BW #cardiac output adjusted for temperature and BW of the simulation (ml/d)
    
    ####################################################################################################
    
    ##Flows to organs (ml/d)
    Fliv    = (liv_frac    * Fcard)
    Fgon    = (gon_frac    * Fcard)
    Fgit    = (git_frac    * Fcard)
    Ffat	  = (fat_frac    * Fcard)
    Fbrain  = (brain_frac  * Fcard)
    Fkidney = (kidney_frac * Fcard)
    Fskin   = (skin_frac   * Fcard)
    Frp     = (rp_frac     * Fcard)
    Fpp     = (Fcard-(Fliv+Fgon+Fgit+Ffat+Fbrain+Fkidney+Fskin+Frp))
    
    Ftotal  = (Fliv+Fgon+Fgit+Ffat+Fbrain+Fkidney+Fskin+Frp+Fpp)
    Fbal    = (Fcard- Ftotal)
    Fegg    = Fgon
    
    VO2     = (VO2_ref*KT(T=TC_k,TR=TR_VO2,TA=TA)*(BW/BW_VO2_ref)^(-0.1))*BW   # O2 consumption rate (mg/d/g) adjusted for T and BW of the simmulation
    Co2w    = ((-0.24 * TC_c + 14.04) * Sat)/(10^3)                            # (mg O2/mL)C of O2 in water at T in celsius
    y_water = VO2/(OEE*Co2w) *( 1/(1000^(0.25) * BW^(0.75)) )  # gill ventilation coefficient OEE=0.71
    
    L0 = A_blood
    
    if(!all(L0==0)){
      PL_old=PL # in M
      K=10.0^(-pKd) # -log10(K) to M
      
      # setting tolerance and maximum iterations
      toler=1E-6
      maxiter = length(K)*1000
      
      L0 = (L0*1E-6)/(MW*(V_blood*1e-3))  # conversion from ng to M
      #PL_old = (PL_old*1E-9)/(MW*V_ven)
      
      PL_new = c(rep(0,length(K)))
      for (k in seq(0,maxiter)){
        for (j in seq(1,length(K))){
          
          # prot available to ligand j
          b = (P_avail + PL_old[j] + L0[j] + K[j])
          c = (P_avail + PL_old[j]) * L0[j]
          #print(glue("b: ", b, "; c: ", c))
          PL_new[j] = (b - sqrt(b * b - 4 * c))/2
          P_avail = P_avail+PL_old[j]-PL_new[j]
        }
        tm = sum(abs((PL_old-PL_new)/PL_new))
        
        if (tm<toler){
          break
        }
        PL_old = PL_new
      }
      PL = PL_new # PL values are stored in M. conversion done in dxdt
    }
    
    blood_bound = PL*MW*1000*V_blood # A_blood_bound in microg
    free = if(all(A_blood)!=0) (1-(blood_bound/A_blood)) else rep(1,n)
    
    dV_urine = urine_rate * BW
    
    ## Compound vectors##
    
    logKow_ion = logKow - delta_Kow
    
    fn_fish = 1/(1+10^(i*(7.4-pKa)))  # (i = 1 for acid and -1 for base) !!! i and pKa - vectors
    
    # fn_gillsurf = 1/(1+10^(i*(6.4-pKa)))
    
    dow = fn_fish * 10^(logKow) + (1-fn_fish)*10^(logKow_ion)   # tmps - vectors
    
    # dow_gillsurf = fn_gillsurf * 10^(logKow) + (1-fn_gillsurf)*10^(logKow_ion)
    
    PBW = 0.008*0.3* dow + 0.007*2.0* dow^(0.94) + 0.134*2.9* dow^(0.63) + 0.851

    y_blood = Fcard * PBW * ( 1/(1000^(0.25) * BW^(0.75)) )   # blood perfusion coefficient vector
    
    kx = (((BW/1000)^(0.75)) /( 2.8*10^(-3) + 68/dow + 1/y_water + 1/y_blood)*1000)
    kout = (1/(68*(dow-1)+1) * kx)

    ### ABSORPTION ###
    
    dA_admin_gill = kx * (A_water/V_water)

    dA_lumen = ( A_bile*K_BG - Ku*A_lumen - Ke_feces*A_lumen)
    
    ### ELIMINATION/ECXRETION ###
    dA_excr_gill = kout *(A_blood* free/V_blood)
    
    
    dA_urine = Ke_urine * A_kidney 
    dA_urine_cum = Ke_urine * A_kidney
    
    dA_feces = Ke_feces * A_lumen
    
    dA_excr_water = ( dA_excr_gill + dA_urine_cum + dA_feces)
    
    dA_bile = Cl_bile * (A_liv/V_liv) - A_bile*K_BG
    
    dA_blood = dA_admin_gill - dA_excr_gill - free*(A_blood/V_blood) * (Fliv + Ffat + Fskin + Fgon + Fgit + Fbrain + Fkidney + Frp + Fpp) +
      free*(Ffat*((A_fat/V_fat)/Pfatb) + a_Fs*Fskin*((A_skin/V_skin)/Pskinb)+ Fbrain*((A_brain/V_brain)/Pbb)+
              (Fliv+Fgon+Fgit+Frp)*((A_liv/V_liv)/Plivb) + (Fkidney+(1-a_Fs)*Fskin+(1-a_Fpp)*Fpp)*((A_kidney/V_kidney)/Pkidb) + a_Fpp*Fpp*((A_pp/V_pp)/Pppb))
    
    dA_liv    = free*(Fliv*(A_blood/V_blood) + Fgon*((A_gon/V_gon)/Pgonb) + Frp*((A_rp/V_rp)/Prpb) + Fgit*((A_git/V_git)/Pgitb)) - 
      free*(Fliv + Frp + Fgit + Fgon)*((A_liv/V_liv)/Plivb)- Cl_bile*(A_liv/V_liv)

    dA_egg    = PS/V_egg* (((A_gon/V_gon)/Pgonb)-((A_egg/V_egg)/Pegggon))
    dA_egg_cum  = PS/V_egg* (((A_gon/V_gon)/Pgonb)-((A_egg/V_egg)/Pegggon))
    
    dA_gon      = free*(Fgon*(A_blood/V_blood) - Fgon*((A_gon/V_gon)/Pgonb)) + PS/V_egg*(((A_egg/V_egg)/Pegggon) - ((A_gon/V_gon)/Pgonb))
    
    dA_fat    = free*Ffat *((A_blood/V_blood) - ((A_fat/V_fat)/Pfatb))
    
    dA_git    = Ku*A_lumen + free*Fgit*((A_blood/V_blood) - (A_git/V_git)/Pgitb)
    
    dA_brain  = free*Fbrain*(A_blood/V_blood - (A_brain/V_brain)/Pbb)
    
    dA_kidney = free*(Fkidney*(A_blood/V_blood) + (1-a_Fpp)*Fpp*((A_pp/V_pp)/Pppb) + (1-a_Fs)*Fskin*((A_skin/V_skin)/Pskinb)) - free*(Fkidney + (1-a_Fpp)*Fpp + (1-a_Fs)*Fskin)*((A_kidney/V_kidney)/Pkidb) - dA_urine
    
    dA_skin   = free*Fskin* (A_blood/V_blood - (A_skin/V_skin)/Pskinb)
    dA_rp     = free*Frp* (A_blood/V_blood - (A_rp/V_rp)/Prpb)
    dA_pp     = free*Fpp* (A_blood/V_blood - (A_pp/V_pp)/Pppb)

    # Chemical in aquarium water
    dA_water = (dA_excr_water - dA_admin_gill)
    
    ### MASS BALANCES ###
    
    A_body = (A_blood + A_bile + A_liv + A_gon + A_brain + A_fat +A_skin + A_kidney + A_git + A_pp + A_rp + A_lumen)
    
    A_elim = (A_excr_water + A_egg_cum)
    
    Mass_bal = (A_admin_gill - A_body - A_elim)
    
    ########## END OF COMPOUNDS #############
    ################### CONCENTRATIONS #############################
    
    ## COMPOUND 1 ##
    C_whole_body = A_body/BW
    C_blood = A_blood/V_blood
    C_plasma = (C_blood)*0.86
    C_liv = A_liv/V_liv
    C_fat = A_fat/V_fat
    C_egg = A_egg/V_egg
    C_gon = A_gon/V_gon
    C_git = A_git/V_git
    C_brain = A_brain/V_brain
    C_kidney = A_kidney/ V_kidney 
    C_skin = A_skin/V_skin
    C_rp  = A_rp/V_rp
    C_pp  = A_pp/V_pp
    C_carcass = (A_body-A_gon-A_liv-A_brain)/(V_total-V_gon-V_liv-V_brain)
    
    ############## END OF CONCENTRATIONS #############
    # 
    list(c(dV_egg,
           dV_urine,
           
           dA_admin_gill, # N3
           dA_lumen,
           dA_excr_gill,
           dA_urine,
           dA_urine_cum,
           dA_feces,
           dA_excr_water,
           dA_bile, # N10
           dA_blood,
           dA_liv,
           dA_egg, #N 13
           dA_egg_cum,
           dA_gon,
           dA_fat,
           dA_git,
           dA_brain,
           dA_kidney,
           dA_skin,
           dA_rp,
           dA_pp,
           dA_water #N23
    ),
    # ,
    "Mass_balance" = Mass_bal,
    "P_avail" = P_avail*1e6,
    "blood_bound" = blood_bound,
    "free" = free,
    "C_whole_body" = C_whole_body,
    "C_blood" = C_blood,
    "C_plasma" = C_plasma,
    "C_liv" = C_liv,
    "C_fat" = C_fat,
    "C_egg" = C_egg,
    "C_gon" = C_gon,
    "C_git" = C_git,
    "C_brain" = C_brain,
    "C_kidney" = C_kidney,
    "C_skin" = C_skin,
    "C_rp" =  C_rp,
    "C_pp" = C_pp,
    "C_carcass" = C_carcass
    )
  }
  )}


#===============================================
# (II.1) Parameters
#===============================================
phys_parms = c(
  # Physiological parameters
  
  Bw_Fcard_ref= 0.5    ,  # Body weight of reference for F_card from Pery, 2014 and Hamilton, 2014
  BW_VO2_ref = 0.4     ,  # Body weight of reference for VO2 from Pery, 2014
  # Environmental condition
  TA           = 3000  ,   # Arrhenius temperature  in Kelvin for zebrafish
  TR_Fcard     = 299.5 , # (Kelvin)
  TR_VO2       = 300.15 ,  # (Kelvin)
  V_water      = 1E12 ,    # Volume of aquarium (mL) large to account for flow-through system
  
  # Effective respiratory volume & cardiac output
  F_card_ref   = 41.328 ,  # reference cardiac output (mL/min/g): Qb_ref = mL/min/g --> 0.0287 * (60*24) = mL/d/g
  VO2_ref      = 9.84, # mg O2/d/g --> calculated from ventilation volume proposed by P?ry, 2014 (Qw= 0.55 mL/min/g)
  OEE          = 0.71, # Oxygen extraction efficiency of 71% proposed by Erickson, 1990
  Sat          = 0.90, # dissolved oxygen saturation of 90% proposed by Erickson, 1990
  art_ven_frac = 1/3, # fraction of arterial blood
  
  sc_blood  = 0.0222, # volume scaling factor : fraction of BW (%)
  sc_gon    = 0.08,
  sc_brain  = 0.01,
  sc_liv    = 0.0205,
  sc_fat    = 0.022 , # fathead minnow data
  sc_skin   = 0.10  , # rainbow trout data
  sc_git    = 0.099 , # rainbow trout data
  sc_kidney = 0.0021,
  sc_rp     = 0.0266,
  sc_pp     = 0.6176  , # (1 -  sc_blood - sc_gon - sc_brain - sc_liv - sc_fat -  sc_skin - sc_git - sc_kidney -sc_rp)
  
  gon_frac    = 0.01252  ,  # Fraction of arterial blood flow
  brain_frac  = 0.0324   ,
  liv_frac    = 0.01942  ,
  fat_frac    = 0.0135   ,
  skin_frac   = 0.05639  ,
  git_frac    = 0.174    ,
  kidney_frac = 0.02233  ,
  rp_frac     = 0.13255  ,
  pp_frac     = 0.53689  , #(1 - gon_frac - brain_frac - liv_frac - fat_frac - skin_frac - git_frac -  kidney_frac - rp_frac)
  
  a_Fpp = 0.4, # Fraction of PPT blood going to venous
  a_Fs  = 0.1, # Fraction of skin blood going to venous
  
  P_avail = 0.0008,
  # Egg #
  Egggrowth	=	0.027	,
  spawnrate = 1.5, # d
  V_one_egg = 0.000212, #V_egg
  delta_Kow = 3.1, 
  
  urine_rate     = 0.057947686, # V_burst = 1.2 mL.kg-1 every 29.82 minutes proposed by Curtis 1991 --> 1.2e-03 mL.g BW-1?
  urination_interval = 30/60/24      # 30min proposed by Curtis 1991 --> 30/60/24hr
)
#----------------------------------------

# Our mixture
exposure_parms <- c(
  TC_c         = 26, # N1
  BW       = 0.8    ,  # Body weight (g fresh weight)
  
  # temps
  start  = 0,       # days
  stop = 32,
  period = NA,       # days between two doses
  time_final_dose = NA,
  time_first_dose = 0,
  time_end_exposure = 32
)
### Chen 
exposure_parms <- c(
  TC_c = 26,
  BW       = 0.4,  # Body weight (g fresh weight) 0.4 - Chen; 0.8 - others
  start  = 0,       # days
  stop = 48,   #changed from 100 simulation length
  period = NA,       # days between two doses
  time_final_dose = NA,
  time_first_dose = 0,
  time_end_exposure = 24
)
### Wen 
exposure_parms <- c(
  TC_c = 26,
  BW       = 0.45,  # Body weight (g fresh weight) 0.4 - Chen; 0.8 - others
  start  = 0,       # days
  stop = 28,   #changed from 100 simulation length
  period = NA,       # days between two doses
  time_final_dose = NA,
  time_first_dose = 0,
  time_end_exposure = 28
)

#-----chemicals with PCs: PFHpA, PFNA, PFHxS, PFOA, PFBS 5 chem-----------------
##### chemicals with PCs: PFHpA, PFNA, PFHxS, PFOA, PFBS 5 chem

chem_parms<-c( # repeated params in this vector: 27
  
  # Chemical parameters
  i=c(1,1,1,1,1),
  pKa = c(0.47, 0.52, -3.34, 0.5,-3.57),  # PFHpA, PFNA, PFHxS, PFOA, PFBS 
  logKow = c(3.94, 4.97, 3.41, 4.35, 3.9), 
  
  MW  = c(364.06,    464.08,    400.12, 414.07,300.1), # PFHpA, PFNA, PFHxS, PFOA, PFBS
  Dmw = c(10^(3.65),10^(3.65),10^(3.65), 10^(3.51), 10^(2.63)), # ADD dummies
  Dpw = c(10^(3.81),10^(3.81),10^(3.81), 10^(3.43), 10^(3.20)), # ADD dummies
  pKd = c(3.36, 3.24, 3.15, 3.08, 2.78), # 
  
  Pegggon = c(0.17, 0.13, 0.19, 0.127, 0.342), # our data 
  Plivb   = c(1.4,  0.26, 0.46, 0.6,   0.9), # average With Chen
  Pgonb   = c(0.13, 0.2,  0.62, 0.4,   0.9), # average With Chen
  Pbb     = c(0.08, 0.07, 0.07, 0.2,   0.09), # average
  Prpb    = c(0.23, 0.25, 0.26, 0.26,  0.27), # Gills 
  Pppb    = c(0.08, 0.05, 0.06, 0.08,  0.09), # muscle
  Pgitb   = c(0.93, 0.31, 0.64, 0.64,  0.7),
  Pfatb   = c(0.5,  0.3,  0.5,  0.5,   0.4),
  Pskinb  = c(0.5,  0.5,  0.5,  0.5,   0.5),
  Pkidb   = c(0.93, 0.31, 0.64, 0.64,  0.7), 
  
  
  #Oral absorption
  Ku       = c(0,0,0,0,0),   # Diffusion coefficient
  frac_abs = c(0,0,0,0,0),
  
  # Excretion
  Cl_bile  = c(0,0,0,0,0), # ml/d
  Ke_urine = c(0,0,0,0,0), # 1/d
  K_BG     = c(0,0,0,0,0), # 1/d
  Ke_feces = c(0.83,0.83,0.83,0.83,0.83) , # 1/d estiamted from Nichols  et al. 2004 0.83
  
  
  PL = c(0,0,0,0,0),
  blood_bound = c(0,0,0,0,0),
  
  # # Exposure quantity (ng/mL = microg/L)

  WaterExposure  = c(9.6*10^(-3), 7.8*10^(-3), 11.4*10^(-3), 8.9*10^(-3), 8.7*10^(-3)) # microg/mL  PFHpA, PFNA, PFHxS, PFOA, PFBS our
)


n=5
yini <- c(V_egg = 0.000212, 
          V_urine = 0, 
          
          A_admin_gill = rep(0,n),
          A_lumen = rep(0,n), 
          A_excr_gill = rep(0,n), 
          A_urine = rep(0,n), 
          A_urine_cum = rep(0,n), 
          A_feces = rep(0,n), 
          A_excr_water = rep(0,n), 
          A_bile = rep(0,n), 
          A_blood = rep(0,n),
          A_liv = rep(0,n),
          A_egg = rep(0,n),
          A_egg_cum = rep(0,n),
          A_gon = rep(0,n), 
          A_fat = rep(0,n), 
          A_git = rep(0,n), 
          A_brain = rep(0,n), 
          A_kidney = rep(0,n), 
          A_skin = rep(0,n), 
          A_rp = rep(0,n), 
          A_pp = rep(0,n),
          A_water = unname(chem_parms)[(1+n*25):((1+n-1)+n*25)]*unname(phys_parms["V_water"])
)
#-------------------------------

phys_exp_parms = c(phys_parms, exposure_parms)

### Events ###
water_exp_event = rep(0,n)
for(i in seq(1,n)){
  if(n>1){
    water_exp_event[i]=glue("A_water", i)}
  if(n==1){
    (water_exp_event[1]="A_water")}
}
spawn_event = rep(0,n)
for(i in seq(1,n)){
  if(n>1){
    spawn_event[i]=glue("A_egg", i)}
  if(n==1){
    (spawn_event[1]="A_egg")}
}

events_end_exposure_w <- list(data = rbind(data.frame(var = water_exp_event, 
                                                      time  = as.numeric(phys_exp_parms["time_end_exposure"]),
                                                      value = 0,
                                                      method = c("replace"))
))

events_spawn <-          list(data = rbind(data.frame(var = rep(spawn_event,each=length(seq(phys_exp_parms["spawnrate"],phys_exp_parms["stop"],phys_exp_parms["spawnrate"]))),
                                                      time  =seq(phys_exp_parms["spawnrate"],phys_exp_parms["stop"],phys_exp_parms["spawnrate"]), 
                                                      # seq(phys_exp_parms["spawnrate"],phys_exp_parms["stop"],phys_exp_parms["spawnrate"])),
                                                      value = 0,
                                                      method = c("replace"))
))

# events<-list(data = rbind(events_end_exposure_w[[1]], events_spawn[[1]]))
events<-list(data = rbind(events_end_exposure_w[[1]]))

# events<-NULL


PFHpA_mix_our = data.frame(time = c(1,3,7,21,32),
                           Carcass = c(12,19.13333333,37.26666667,35.35,37.33333333)/1000,
                           Ovaries = c(0.5,17.03333333,32.94,21.72,42.7)/1000,
                           sd_carcass = c(1.587450787,3.219213154,7.890711839,6.151828996,12.97009381)/1000,
                           sd_ovaries = c(0,4.933896364,23.73008217,14.133761,13.5114766)/1000) # microg/g

PFHxS_mix_our = data.frame(time = c(1,3,7,21,32),
                           Carcass = c(26.83333333,46.23333333,120.6666667,186,224.6666667)/1000,
                           Ovaries = c(43.9,102.4666667,314,708.6666667,1122)/1000,
                           Liver = c(81.73333333,145.3333333,333.6666667,638.3333333,685)/1000,
                           sd_carcass = c(3.855299383,4.245389656,23.71356855,29.13760457,47.50087718)/1000,
                           sd_ovaries = c(3.764306045,13.45263295,51.97114584,81.11925378,228.6919325)/1000,
                           sd_liver = c(23.00028985,9.291573243,30.66485502,207.403793,118.8822947)/1000) # microg/g

PFOA_mix_our = data.frame(time = c(1,3,7,21,32),
                          Carcass = c(37.33, 67.03, 177.00, 259.67,292.33)/1000,
                          Ovaries = c(23.07, 52.17, 139.00, 314.67, 479.00)/1000,
                          Brain = c(2.69, 31.80, 112.00, 171.67, 412.67)/1000,
                          Liver = c(38.83, 95.63, 208.00, 400.33, 453.67)/1000,
                          sd_carcass = c(4.54, 17.17, 59.23, 101.50, 80.05)/1000,
                          sd_ovaries = c(1.72, 11.36, 38.69, 78.49, 122.75)/1000,
                          sd_brain = c(3.79, 21.08, 49.33, 111.14, 285.56)/1000,
                          sd_liver = c(31.91, 26.39, 56.29, 87.08,166.50)/1000) # microg/g

PFNA_mix_our = data.frame(time = c(1,3,7,21,32),
                          Carcass = c(121.5666667,250,608.3333333,1275,1754)/1000,
                          Ovaries = c(90.6,157.3333333,666.6666667,2619,3617.666667)/1000,
                          Brain = c(120.5333333,208.3333333,494.6666667,1042.666667,1099.666667)/1000,
                          Liver = c(127,259.6666667,813,2327,2242)/1000,
                          sd_carcass = c(33.45688768,42.53234064,56.03867712,130.8472392,245.2203091)/1000,
                          sd_ovaries = c(18.45020325,27.68272626,59.18051481,356.4210993,637.1729226)/1000,
                          sd_brain = c(63.32182352,26.1023626,108.9602374,251.8518083,85.40686936)/1000,
                          sd_liver = c(51.39066063,26.38812864,159.5587666,873.0515449,545.0165135)/1000) # microg


PFBS_mix_our = data.frame(time = c(1,3,7,21,32),
                          Carcass = c(5.563333333,17.23,16.6,12.23333333,10.45)/1000,
                          Ovaries = c(13.83333333,34.46666667,44.83333333,43.83333333,51.7)/1000,
                          sd_carcass = c(0.927918818,12.1663347,6.921704992,1.171893055,2.795371174)/1000,
                          sd_ovaries = c(1.327905619,2.218858565,22.80182741,21.70860045,18.70614872)/1000) # microg/g

Times = sort(unique(c(seq(phys_exp_parms["start"],phys_exp_parms["stop"]), PFOA_mix_our$time)))

#===============================================
# (III) Numerical resolution and post-calcul
#===============================================

out <- ode(times = Times,
           func  = ZF.model.c,
           y     = yini ,
           parms = phys_exp_parms,
           chem_parms = chem_parms,
           method="lsodes",
           events = events)

