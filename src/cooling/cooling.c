/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/cooling/cooling.c
 * \date        05/2018
 * \brief       Module for gas radiative cooling
 * \details     contains functions:
 *                double DoCooling(double u_old, double rho, double dt, double
 *                  *ne_guess)
 *                double GetCoolingTime(double u_old, double rho, double
 *                  *ne_guess)
 *                double convert_u_to_temp(double u, double rho, double
 *                  *ne_guess)
 *                void find_abundances_and_rates(double logT, double rho,
 *                  double *ne_guess)
 *                double CoolingRateFromU(double u, double rho, double
 *                  *ne_guess)
 *                void SetOutputGasState(int i, double *ne_guess, double *nH0,
 *                  double *coolrate)
 *                double CoolingRate(double logT, double rho, double *nelec)
 *                void MakeRateTable(void)
 *                void ReadIonizeParams(char *fname, int which)
 *                void IonizeParamsUVB(void)
 *                void SetZeroIonization(void)
 *                void IonizeParams(void)
 *                void InitCool(void)
 *                void cooling_only(void)
 *                void cool_cell(int i)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

 #include <math.h>
 #include <mpi.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <hdf5.h>
 #include "../main/allvars.h"
 #include "../main/proto.h"
 
 #ifdef COOLING
 
 static double Tmin = 0.0;     /*!< min temperature in log10 */
 static double Tmax = 9.0;     /*!< max temperature in log10 */
 static double deltaT;         /*!< log10 of temperature spacing in the interpolation tables */
 static GasState gs;           /*!< gas state */
 static RateTable *RateT;      /*!< tabulated rates */
 static PhotoTable *PhotoTUVB; /*!< photo-ionization/heating rate table for UV background */
 static PhotoCurrent pc;       /*!< current interpolated photo rates */
 static int NheattabUVB;       /*!< length of UVB photo table */
 static DoCoolData DoCool;     /*!< cooling data */
#ifdef METALLIC_COOLING
 static MetalTable MetalT;  /*!< heating functions and values the high temperature metallic cooling */
 static MetalTable MetalT_CIE;  /*!< heating functions and values the high temperature metallic cooling */

#endif 
 /*! \brief Computes the new internal energy per unit mass.
  *
  *  The function solves for the new internal energy per unit mass of the gas
  *  by integrating the equation for the internal energy with an implicit
  *  Euler scheme. The root of resulting non linear equation,
  *  which gives tnew internal energy, is found with the bisection method.
  *  Arguments are passed in code units.
  *
  *  \param[in] u_old the initial (before cooling is applied) internal energy
  *             per unit mass of the gas cell.
  *  \param[in] rho   the proper density of the gas cell.
  *  \param[in] dt    the duration of the time step.
  *  \param[in] ne_guess electron number density relative to hydrogen number
  *             density (for molecular weight computation).
  *
  *  \return The new internal energy per unit mass of the gas cell.
  */
 double DoCooling(double u_old, double rho, double dt, double *ne_guess, double metallicity) 
 {
   double u, du;
   double u_lower, u_upper;
   double ratefact;
   double LambdaNet;
 
   int iter = 0;
 
   DoCool.u_old_input    = u_old;
   DoCool.rho_input      = rho;
   DoCool.dt_input       = dt;
   DoCool.ne_guess_input = *ne_guess;
 
   if(!gsl_finite(u_old))
     terminate("invalid input: u_old=%g\n", u_old);
 
   if(u_old < 0 || rho < 0)
     terminate("invalid input: task=%d u_old=%g  rho=%g  dt=%g  All.MinEgySpec=%g\n", ThisTask, u_old, rho, dt, All.MinEgySpec);
 
   rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
   u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
   dt *= All.UnitTime_in_s / All.HubbleParam;
 
   gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
   ratefact = gs.nHcgs * gs.nHcgs / rho;
 
   u       = u_old;
   u_lower = u;
   u_upper = u;
 
   LambdaNet = CoolingRateFromU(u, rho, ne_guess, metallicity);
 
   /* bracketing */
   if(u - u_old - ratefact * LambdaNet * dt < 0) /* heating */
     {
       u_upper *= sqrt(1.1);
       u_lower /= sqrt(1.1);
       while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess, metallicity) * dt < 0)
         {
           u_upper *= 1.1;
           u_lower *= 1.1;
         }
     }
 
   if(u - u_old - ratefact * LambdaNet * dt > 0)
     {
       u_lower /= sqrt(1.1);
       u_upper *= sqrt(1.1);
       while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess, metallicity) * dt > 0)
         {
           u_upper /= 1.1;
           u_lower /= 1.1;
         }
     }
 
   do
     {
       u = 0.5 * (u_lower + u_upper);
 
       LambdaNet = CoolingRateFromU(u, rho, ne_guess, metallicity);
 
       if(u - u_old - ratefact * LambdaNet * dt > 0)
         {
           u_upper = u;
         }
       else
         {
           u_lower = u;
         }
 
       du = u_upper - u_lower;
 
       iter++;
 
       if(iter >= (MAXITER - 10))
         printf("u= %g\n", u);
     }
   while(fabs(du / u) > 1.0e-6 && iter < MAXITER);
 
   if(iter >= MAXITER)
     terminate(
         "failed to converge in DoCooling(): DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= "
         "%g\n",
         DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
 
   u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; /* to internal units */
 
   return u;
 }
 
 /*! \brief Returns the cooling time.
  *
  *  If we actually have heating, a cooling time of 0 is returned.
  *
  *  \param[in] u_old The initial (before cooling is applied) internal energy
  *             per unit mass of the gas cell.
  *  \param[in] rho The proper density of the gas cell.
  *  \param[in] ne_guess Electron number density relative to hydrogen number
  *             density (for molecular weight computation).
  *
  *  \return Cooling time; 0 if heating.
  */
 double GetCoolingTime(double u_old, double rho, double *ne_guess, double metallicity)
 {
   double u;
   double ratefact;
   double LambdaNet, coolingtime;
 
   DoCool.u_old_input    = u_old;
   DoCool.rho_input      = rho;
   DoCool.ne_guess_input = *ne_guess;
 
   rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
   u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
 
   gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
   ratefact = gs.nHcgs * gs.nHcgs / rho;
 
   u = u_old;
 
   LambdaNet = CoolingRateFromU(u, rho, ne_guess, metallicity);
 
   if(LambdaNet >= 0) /* ups, we have actually heating due to UV background */
     return 0;
 
   coolingtime = u_old / (-ratefact * LambdaNet);
 
   coolingtime *= All.HubbleParam / All.UnitTime_in_s;
 
   return coolingtime;
 }
 
 /*! \brief Compute gas temperature from internal energy per unit mass.
  *
  *   This function determines the electron fraction, and hence the mean
  *   molecular weight. With it arrives at a self-consistent temperature.
  *   Element abundances and the rates for the emission are also computed.
  *
  *  \param[in] u   internal energy per unit mass.
  *  \param[in] rho gas density.
  *  \param[in, out] ne_guess electron number density relative to hydrogen
  *                  number density
  *
  *  \return The gas temperature.
  */
 double convert_u_to_temp(double u, double rho, double *ne_guess)
 {
   double temp, temp_old, temp_new, max = 0, ne_old;
   double mu;
   int iter = 0;
 
   double u_input, rho_input, ne_input;
 
   u_input   = u;
   rho_input = rho;
   ne_input  = *ne_guess;
 
   mu   = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess);
   temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
  //  if (temp == 0)
  //    terminate("T=0 error in convert_u_to_temp() before we find the abundances, T = %g, mu = %g, u_input = %g\n", temp, mu, u_input);
 
   do
     {
       ne_old = *ne_guess;
       find_abundances_and_rates(log10(temp), rho, ne_guess);
       temp_old = temp;
       
       mu = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess);
 
       temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
      //  if (temp_new == 0)
      //  {
      //    printf("temp_new=%g\n", temp_new);
      //    terminate("temp_new = 0, in convert_u_to_temp: T = %g, u= %g, mu=%g\n", temp, u, mu);
      //  }
       max = dmax(max, temp_new / (1 + gs.yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));
 
       temp = temp_old + (temp_new - temp_old) / (1 + max);
       iter++;
 
       if(iter > (MAXITER - 10))
         printf("-> temp= %g ne=%g\n", temp, *ne_guess);
     }
   while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);
 
   if(iter >= MAXITER)
     {
       printf("failed to converge in convert_u_to_temp()\n");
       printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
       printf("DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n", DoCool.u_old_input,
              DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
       terminate("convergence failure");
     }
 
   gs.mu = mu;
 
   return temp;
 }
 
 /*! \brief Computes the actual abundance ratios.
  *
  *  The chemical composition of the gas is primordial (no metals are present).
  *
  *  \param[in] logT log10 of gas temperature.
  *  \param[in] rho Gas density.
  *  \param[in, out] ne_guess Electron number density relative to hydrogen
  *                  number density.
  *
  *  \return void
  */
 void find_abundances_and_rates(double logT, double rho, double *ne_guess)
 {
   double neold, nenew;
   int j, niter;
   double flow, fhi, t;
 
   double logT_input, rho_input, ne_input;
 
   logT_input = logT;
   rho_input  = rho;
   ne_input   = *ne_guess;
 
  //  if(!gsl_finite(logT))
  //  {
  //    // Per Illustris T = (gamma - 1)*u/kb * UnitMass/UnitEnergy *mu. The only candidates are: u and mu since the rest are constant
  //    printf("Temperature T = %g\n", pow(10.0, logT));
  //    terminate("logT=%g\n", logT);
  //  }
   if(logT <= Tmin) /* everything neutral */
     {
       gs.nH0    = 1.0;
       gs.nHe0   = gs.yhelium;
       gs.nHp    = 0;
       gs.nHep   = 0;
       gs.nHepp  = 0;
       gs.ne     = 0;
       *ne_guess = 0;
       return;
     }
 
   if(logT >= Tmax) /* everything is ionized */
     {
       gs.nH0    = 0;
       gs.nHe0   = 0;
       gs.nHp    = 1.0;
       gs.nHep   = 0;
       gs.nHepp  = gs.yhelium;
       gs.ne     = gs.nHp + 2.0 * gs.nHepp;
       *ne_guess = gs.ne; /* note: in units of the hydrogen number density */
       return;
     }
 
   t    = (logT - Tmin) / deltaT;
   j    = (int)t;
   fhi  = t - j;
   flow = 1 - fhi;
 
   if(*ne_guess == 0)
     *ne_guess = 1.0;
 
   gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
 
   gs.ne    = *ne_guess;
   neold    = gs.ne;
   niter    = 0;
   gs.necgs = gs.ne * gs.nHcgs;
 
   /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
   do
     {
       niter++;
 
       gs.aHp   = flow * RateT[j].AlphaHp + fhi * RateT[j + 1].AlphaHp;
       gs.aHep  = flow * RateT[j].AlphaHep + fhi * RateT[j + 1].AlphaHep;
       gs.aHepp = flow * RateT[j].AlphaHepp + fhi * RateT[j + 1].AlphaHepp;
       gs.ad    = flow * RateT[j].Alphad + fhi * RateT[j + 1].Alphad;
       gs.geH0  = flow * RateT[j].GammaeH0 + fhi * RateT[j + 1].GammaeH0;
       gs.geHe0 = flow * RateT[j].GammaeHe0 + fhi * RateT[j + 1].GammaeHe0;
       gs.geHep = flow * RateT[j].GammaeHep + fhi * RateT[j + 1].GammaeHep;
 
       if(gs.necgs <= 1.e-25 || pc.J_UV == 0)
         {
           gs.gJH0ne = gs.gJHe0ne = gs.gJHepne = 0;
         }
       else
         {
           gs.gJH0ne  = pc.gJH0 / gs.necgs;
           gs.gJHe0ne = pc.gJHe0 / gs.necgs;
           gs.gJHepne = pc.gJHep / gs.necgs;
         }
 
       gs.nH0 = gs.aHp / (gs.aHp + gs.geH0 + gs.gJH0ne); /* eqn (33) */
       gs.nHp = 1.0 - gs.nH0;                            /* eqn (34) */
 
       if((gs.gJHe0ne + gs.geHe0) <= SMALLNUM) /* no ionization at all */
         {
           gs.nHep  = 0.0;
           gs.nHepp = 0.0;
           gs.nHe0  = gs.yhelium;
         }
       else
         {
           gs.nHep =
               gs.yhelium / (1.0 + (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne) + (gs.geHep + gs.gJHepne) / gs.aHepp); /* eqn (35) */
           gs.nHe0  = gs.nHep * (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne);                                          /* eqn (36) */
           gs.nHepp = gs.nHep * (gs.geHep + gs.gJHepne) / gs.aHepp;                                                   /* eqn (37) */
         }
 
       neold = gs.ne;
 
       gs.ne    = gs.nHp + gs.nHep + 2 * gs.nHepp; /* eqn (38) */
       gs.necgs = gs.ne * gs.nHcgs;
 
       if(pc.J_UV == 0)
         break;
 
       nenew    = 0.5 * (gs.ne + neold);
       gs.ne    = nenew;
       gs.necgs = gs.ne * gs.nHcgs;
 
       if(fabs(gs.ne - neold) < 1.0e-4)
         break;
 
       if(niter > (MAXITER - 10))
         printf("ne= %g  niter=%d\n", gs.ne, niter);
     }
   while(niter < MAXITER);
 
   if(niter >= MAXITER)
     {
       printf("gs.aHp = %le\n", gs.aHp);
       char buff[1000];
       sprintf(buff, "%s/cooling_task%d.dat", All.OutputDir, ThisTask);
       FILE *fp = fopen(buff, "w");
       fwrite(&All.Time, sizeof(double), 1, fp);
       fwrite(&logT_input, sizeof(double), 1, fp);
       fwrite(&rho_input, sizeof(double), 1, fp);
       fwrite(&ne_input, sizeof(double), 1, fp);
       fclose(fp);
       terminate(
           "no convergence reached in find_abundances_and_rates(): logT_input= %g  rho_input= %g  ne_input= %g "
           "DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n",
           logT_input, rho_input, ne_input, DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
     }
   gs.bH0  = flow * RateT[j].BetaH0 + fhi * RateT[j + 1].BetaH0;
   gs.bHep = flow * RateT[j].BetaHep + fhi * RateT[j + 1].BetaHep;
   gs.bff  = flow * RateT[j].Betaff + fhi * RateT[j + 1].Betaff;
 
   *ne_guess = gs.ne;
 }
 
 /*! \brief Get cooling rate from gas internal energy.
  *
  *  This function first computes the self-consistent temperature
  *  and abundance ratios, and then it calculates
  *  (heating rate-cooling rate)/n_h^2 in cgs units.
  *
  *  \param[in] u Gas internal energy per unit mass.
  *  \param[in] rho Gas density.
  *  \param[in, out] ne_guess Electron number density relative to hydrogen
  *                  number density.
  *
  *  \return Cooling rate.
  */
 double CoolingRateFromU(double u, double rho, double *ne_guess, double metallicity)
 {
   double temp;
  //  if (u == 0)
  //  {
  //   terminate("In CoolingRateFromU, u = %g\n", u); 
  //  }
   temp = convert_u_to_temp(u, rho, ne_guess);
 
   return CoolingRate(log10(temp), rho, ne_guess, metallicity);
 }
 
 /*! \brief  This function computes the self-consistent temperature and
  *          abundance ratios.
  *
  *  Used only in io_fields.c for calculating output fields.
  *
  *  \param[in] i index into SphP for gas cell to consider.
  *  \param[in, out] ne_guess pointer to electron number density relative to
  *                  hydrogen number density (modified).
  *  \param[out] nH0 Pointer to the neutral hydrogen fraction (set to value in
  *              the GasState struct).
  *  \param[out] coolrate Pointer to cooling rate (set to value from
  *              CoolingRateFromU).
  *
  *  \return void
  */
 void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate, double *z_lambda, double metallicity)
 {
   double sfr = 0;
   double rho = SphP[i].Density * All.cf_a3inv;
   double u   = dmax(All.MinEgySpec, SphP[i].Utherm);
 
   /* update GasState as appropriate given compile-time options and cell properties */
   //  #if defined(USE_SFR)
   //  sfr = get_starformation_rate(i);  // call is superfluous at this place
   // #endif
 
   /* update DoCool */
   DoCool.u_old_input    = u;
   DoCool.rho_input      = rho;
   DoCool.ne_guess_input = *ne_guess;
 
   /* convert to physical cgs units */
   rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
   u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
 
   /* calculate cooling rate (and so ne_guess and all of gs including nH0, nHeII) */
   *coolrate = CoolingRateFromU(u, rho, ne_guess, metallicity);
#ifdef METALLIC_COOLING
   double T = convert_u_to_temp(u, rho, ne_guess);
#ifdef CIE_PIE_COOLING
  *z_lambda = LambdaMetals(T, gs.nHcgs, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
#endif 
#ifdef CIE_COOLING
   *z_lambda = LambdaMetals_CIE(T, MetalT_CIE.lambda_m, MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
#endif
#endif /* #ifdef METALLIC_COOLING */
   *nH0 = gs.nH0;
 }
 
 /*! \brief  Calculate (heating rate-cooling rate)/n_h^2 in cgs units.
  *
  *  \param[in] logT log10 of gas temperature.
  *  \param[in] rho Gas density.
  *  \param[in, out] nelec Electron number density relative to hydrogen number
  *                  density.
  *
  *  \return (heating rate-cooling rate)/n_h^2.
  */
 double CoolingRate(double logT, double rho, double *nelec, double metallicity)
 {
   double Lambda, Heat;
   double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
   double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
   double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
   double redshift;
   double T;
   double LambdaPrim = 0.0, LambdaMet = 0.0, LambdaDust = 0.0, LambdaMol = 0.0;
 
   if(logT <= Tmin)
     logT = Tmin + 0.5 * deltaT; /* floor at Tmin */
     
    //  if(!gsl_finite(logT))
    //  {
    //    printf("logT=%g\n", logT);
    //    terminate("High T in CoolingRate() below logT <= Tmin, T = %g\n", T); // NOTE: This should be 0 if I understand this correctly.
    //  }
 
   gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
 
   if(logT < Tmax)
     {
 
       find_abundances_and_rates(logT, rho, nelec);
 
       /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
       
       // if(!gsl_finite(u_old))
       // terminate("invalid input: u_old=%g\n", u_old);
 
       T = pow(10.0, logT);
      //  if (T == 0)
      //  {
      //    printf("logT=%g\n", logT);
      //    terminate("T == 0 in CoolingRate() below /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */. T = %g\n", T);
      //  }
       LambdaExcH0   = gs.bH0 * gs.ne * gs.nH0;
       LambdaExcHep  = gs.bHep * gs.ne * gs.nHep;
       LambdaExc     = LambdaExcH0 + LambdaExcHep; /* excitation */
       LambdaIonH0   = 2.18e-11 * gs.geH0 * gs.ne * gs.nH0;
       LambdaIonHe0  = 3.94e-11 * gs.geHe0 * gs.ne * gs.nHe0;
       LambdaIonHep  = 8.72e-11 * gs.geHep * gs.ne * gs.nHep;
       LambdaIon     = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep; /* ionization */
       LambdaRecHp   = 1.036e-16 * T * gs.ne * (gs.aHp * gs.nHp);
       LambdaRecHep  = 1.036e-16 * T * gs.ne * (gs.aHep * gs.nHep);
       LambdaRecHepp = 1.036e-16 * T * gs.ne * (gs.aHepp * gs.nHepp);
       LambdaRecHepd = 6.526e-11 * gs.ad * gs.ne * gs.nHep;
       LambdaRec     = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;
       LambdaFF      = gs.bff * (gs.nHp + gs.nHep + 4 * gs.nHepp) * gs.ne;
       LambdaPrim    = LambdaExc + LambdaIon + LambdaRec + LambdaFF;
 
       if(All.ComovingIntegrationOn)
         {
           redshift    = 1 / All.Time - 1;
           LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
         }
       else
         LambdaCmptn = 0;
         
#ifdef METALLIC_COOLING
      // takes temperature in non-logged T, 
#ifdef CIE_PIE_COOLING
      double lambda_metals = LambdaMetals(T, gs.nHcgs, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
      double electron_solar =  neSolar(T, gs.nHcgs,  MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
#endif
#ifdef CIE_COOLING
      double lambda_metals =  LambdaMetals_CIE(T, MetalT_CIE.lambda_m, MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
      double electron_solar =  neSolar_CIE(T,  MetalT.ne_solar, MetalT.temperature_bins, MetalT.rows_cooling);
#endif

      if (electron_solar <= 0)
        terminate("The solar electron abundance is less than or equal to 0 for: T = %g, nh=%g, ne_solar=%g\n", T, gs.nHcgs, electron_solar);

      double divisor = gs.ne/electron_solar;

      LambdaMet += lambda_metals*metallicity/gs.XZ*divisor; //electron_solar; // gs.ne/electron_solar*

#endif 
         
       Lambda = LambdaPrim + LambdaMet + LambdaDust + LambdaCmptn + LambdaMol;
 
       Heat = 0;
       if(pc.J_UV != 0)
         Heat += (gs.nH0 * pc.epsH0 + gs.nHe0 * pc.epsHe0 + gs.nHep * pc.epsHep) / gs.nHcgs;
     }
   else /* here we're outside of tabulated rates, T>Tmax K */
     {
       /* at high T (fully ionized); only free-free and Compton cooling are present. Assumes no heating. */
       Heat = 0;
 
       LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp =
           LambdaRecHepd                                                                                   = 0;
 
       /* very hot: H and He both fully ionized */
       gs.nHp   = 1.0;
       gs.nHep  = 0;
       gs.nHepp = gs.yhelium;
       gs.ne    = gs.nHp + 2.0 * gs.nHepp;
       *nelec   = gs.ne; /* note: in units of the hydrogen number density */
 
       T        = pow(10.0, logT);
      //  if(!gsl_finite(T))
      //  {
      //    printf("logT=%g\n", logT);
      //    terminate("High T in CoolingRate() where we're outside the tabulated rates: T = %g\n", T);
      //  }
       LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (gs.nHp + 4 * gs.nHepp) * gs.ne;
 
       if(All.ComovingIntegrationOn)
         {
           redshift = 1 / All.Time - 1;
           /* add inverse Compton cooling off the microwave background */
           LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
         }
       else
         LambdaCmptn = 0;
 
       Lambda = LambdaFF + LambdaCmptn;
     }
#if defined(CIE_COOLING)
   return (- Lambda); // CIE ignores the effects of the UV background. So there is no photoionic heating
#endif 
   return (Heat - Lambda);
 }
 
 /*! \brief Make cooling rates interpolation table.
  *
  *  Set up interpolation tables in T for cooling rates given in
  *  KWH, ApJS, 105, 19.
  *
  *  \return void
  */
 void MakeRateTable(void)
 {
   int i;
   double T;
   double Tfact;
 
   gs.yhelium = (1 - gs.XH) / (4 * gs.XH);
   gs.mhboltz = PROTONMASS / BOLTZMANN;
   if(All.MinGasTemp > 0.0)
     Tmin = log10(0.1 * All.MinGasTemp);
   else
     Tmin = 1.0;
   deltaT    = (Tmax - Tmin) / NCOOLTAB;
   gs.ethmin = pow(10.0, Tmin) * (1. + gs.yhelium) / ((1. + 4. * gs.yhelium) * gs.mhboltz * GAMMA_MINUS1);
   /* minimum internal energy for neutral gas */
 
   for(i = 0; i <= NCOOLTAB; i++)
     {
       RateT[i].BetaH0 = RateT[i].BetaHep = RateT[i].Betaff = RateT[i].AlphaHp = RateT[i].AlphaHep = RateT[i].AlphaHepp =
           RateT[i].Alphad = RateT[i].GammaeH0 = RateT[i].GammaeHe0 = RateT[i].GammaeHep = 0;
 
       T     = pow(10.0, Tmin + deltaT * i);
       Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
 
       /* collisional excitation */
       /* Cen 1992 */
       if(118348 / T < 70)
         RateT[i].BetaH0 = 7.5e-19 * exp(-118348 / T) * Tfact;
       if(473638 / T < 70)
         RateT[i].BetaHep = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;
 
       /* free-free */
       RateT[i].Betaff = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
 
       /* recombination */
       /* Cen 1992 */
       /* Hydrogen II */
       RateT[i].AlphaHp = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
       /* Helium II */
       RateT[i].AlphaHep = 1.5e-10 * pow(T, -0.6353);
       /* Helium III */
       RateT[i].AlphaHepp = 4. * RateT[i].AlphaHp;
 
       /* Cen 1992 */
       /* dielectric recombination */
       if(470000 / T < 70)
         RateT[i].Alphad = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));
 
       /* collisional ionization */
       /* Cen 1992 */
       /* Hydrogen */
       if(157809.1 / T < 70)
         RateT[i].GammaeH0 = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
       /* Helium */
       if(285335.4 / T < 70)
         RateT[i].GammaeHe0 = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
       /* Hellium II */
       if(631515.0 / T < 70)
         RateT[i].GammaeHep = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
     }
 }
 
 /*! \brief Read table input for ionizing parameters.
  *
  *  \param[in] fname Name of file that contains the tabulated parameters.
  *  \param[in] which Flag used to identify the type of the ionizing background
  *                   (0 = UV background, 1 = AGN background, 2=RADCOOL).
  *
  *  \return void
  */
 void ReadIonizeParams(char *fname, int which)
 {
   int iter, i;
   FILE *fdcool;
   float dummy;
 
   if(which == 0)
     {
       NheattabUVB = 0;
 
       for(iter = 0, i = 0; iter < 2; iter++)
         {
           if(!(fdcool = fopen(fname, "r")))
             terminate("COOLING: cannot read ionization table in file `%s'\n", fname);
           if(iter == 0)
             while(fscanf(fdcool, "%g %g %g %g %g %g %g", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
               NheattabUVB++;
           if(iter == 1)
             while(fscanf(fdcool, "%g %g %g %g %g %g %g", &PhotoTUVB[i].variable, &PhotoTUVB[i].gH0, &PhotoTUVB[i].gHe,
                          &PhotoTUVB[i].gHep, &PhotoTUVB[i].eH0, &PhotoTUVB[i].eHe, &PhotoTUVB[i].eHep) != EOF)
               i++;
           fclose(fdcool);
 
           if(iter == 0)
             {
               PhotoTUVB = (PhotoTable *)mymalloc("PhotoT", NheattabUVB * sizeof(PhotoTable));
               mpi_printf("COOLING: read ionization table with %d entries in file `%s'.\n", NheattabUVB, fname);
             }
         }
       /* ignore zeros at end of treecool file */
       for(i = 0; i < NheattabUVB; ++i)
         if(PhotoTUVB[i].gH0 == 0.0)
           break;
 
       NheattabUVB = i;
       mpi_printf("COOLING: using %d ionization table entries from file `%s'.\n", NheattabUVB, fname);
     }
 }
 
 /*! \brief Set the ionization parameters for the UV background.
  *
  *  \return void
  */
 void IonizeParamsUVB(void)
 {
   int i, ilow;
   double logz, dzlow, dzhi;
   double redshift;
 
   if(All.ComovingIntegrationOn)
     redshift = 1 / All.Time - 1;
   else
     {
       redshift = 0.0;
     }
 
   logz = log10(redshift + 1.0);
   ilow = 0;
   for(i = 0; i < NheattabUVB; i++)
     {
       if(PhotoTUVB[i].variable < logz)
         ilow = i;
       else
         break;
     }
 
   dzlow = logz - PhotoTUVB[ilow].variable;
   dzhi  = PhotoTUVB[ilow + 1].variable - logz;
 
   if(NheattabUVB == 0 || logz > PhotoTUVB[NheattabUVB - 1].variable || PhotoTUVB[ilow].gH0 == 0 || PhotoTUVB[ilow + 1].gH0 == 0)
     {
       SetZeroIonization();
       return;
     }
   else
     pc.J_UV = 1;
 
   pc.gJH0   = pow(10., (dzhi * log10(PhotoTUVB[ilow].gH0) + dzlow * log10(PhotoTUVB[ilow + 1].gH0)) / (dzlow + dzhi));
   pc.gJHe0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHe) + dzlow * log10(PhotoTUVB[ilow + 1].gHe)) / (dzlow + dzhi));
   pc.gJHep  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHep) + dzlow * log10(PhotoTUVB[ilow + 1].gHep)) / (dzlow + dzhi));
   pc.epsH0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].eH0) + dzlow * log10(PhotoTUVB[ilow + 1].eH0)) / (dzlow + dzhi));
   pc.epsHe0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHe) + dzlow * log10(PhotoTUVB[ilow + 1].eHe)) / (dzlow + dzhi));
   pc.epsHep = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHep) + dzlow * log10(PhotoTUVB[ilow + 1].eHep)) / (dzlow + dzhi));
 
   return;
 }
 
 /*! \brief Reset the ionization parameters.
  *
  *  \return void
  */
 void SetZeroIonization(void) { memset(&pc, 0, sizeof(PhotoCurrent)); }
 
 /*! \brief Wrapper function to set the ionizing background.
  *
  *  \return void
  */
 void IonizeParams(void) { IonizeParamsUVB(); }
 
 /*! \brief Initialize the cooling module.
  *
  *  This function initializes the cooling module. In particular,
  *  it allocates the memory for the cooling rate and ionization tables
  *  and initializes them.
  *
  *  \return void
  */
 void InitCool(void)
 {
   /* set default hydrogen mass fraction */
   gs.XH = HYDROGEN_MASSFRAC;
 
   /* zero photo-ionization/heating rates */
   SetZeroIonization();
 
   /* allocate and construct rate table */
   RateT = (RateTable *)mymalloc("RateT", (NCOOLTAB + 1) * sizeof(RateTable));
   ;
   MakeRateTable();
 
   /* read photo tables */
   ReadIonizeParams(All.TreecoolFile, 0);
 
#ifdef METALLIC_COOLING
   /* Read the cooling table */
  // Taken from (Wiersma, R. P. C., Schaye, J., & Smith, B. D. 2009a, MNRAS, 393, 99) 
  // The hdf5 file contains cooling rates Lambda_{odot}/nh^2 for various metals for given temperature and number density bins and for solar abundances
   ReadMetallicParams("./z_0.000.hdf5"); 
   ReadMetallicParams_CIE("./z_collis.hdf5"); 

   /* Set default metalicity mass fraction */
   // The values here are taken from (Hopkins 2018, MNRAS, 480, 800) and are taken for solar abundances
   // (Z, He, C, N, O, Ne, Mg, Si, S, Ca, Fe)
   // (0.02, 0.28, 3.26e-3, 1.3e-3, 8.65e-3, 2.22e-3, 9.31e-4, 1.08e-3, 6.44e-4, 1.01e-4, 1.73e-3)
   gs.XZ = 0.02; // gs.XC = 3.26e-3, gs.XN=1.3e-3, gs.XO=8.65e-3, gs.XNe=2.22e-3, gs.XMg=9.31e-4, gs.XSi=1.08e-3, gs.XS=6.44e-4, gs.XCa=1.01e-4, gs.XFe=1.73e-3;
   
  //  double T1 = 1e5;
  //  double T2 = 1.874e06;
  //  double T3 = 5e5;
  //  double T4 = -1;
  //  double T5 = 0;
  //  double T6 = 1e30;
  //  double T7 = 1e4;

  //   // //  testing functions
  //  double lambda_1 = LambdaMetals_CIE(T1,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T1, lambda_1);
  //  double lambda_2 = LambdaMetals_CIE(T2,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T2, lambda_2);
  //  double lambda_3 = LambdaMetals_CIE(T3,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T3, lambda_3);
  //  double lambda_4 = LambdaMetals_CIE(T4,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T4, lambda_4);
  //  double lambda_5 = LambdaMetals_CIE(T5,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T5, lambda_5);
  //  double lambda_6 = LambdaMetals_CIE(T6,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T6, lambda_6);
  //  double lambda_7 = LambdaMetals_CIE(T7,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T7, lambda_7);
  //  double lambda_8 = LambdaMetals_CIE(T7,MetalT_CIE.lambda_m,MetalT_CIE.temperature_bins, MetalT_CIE.rows_cooling);
  //  printf("Interpolated cooling for T = %0.3e, Lambda_z = %0.3e\n", T7, lambda_7);
   

  //  printf("testing functions\n");
  // //  testing functions
  //  double lambda_1 = LambdaMetals(min_T, max_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", min_T, max_rho, lambda_1);
  //  double lambda_2 = LambdaMetals(max_T, min_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", max_T, min_rho, lambda_2);
  //  double lambda_3 = LambdaMetals(min_T, min_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", min_T, min_rho, lambda_3);
  //  double lambda_4 = LambdaMetals(max_T, max_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", max_T, max_rho, lambda_4);
  //  double lambda_5 = LambdaMetals(low_T, low_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", low_T, low_rho, lambda_5);
  //  double lambda_6 = LambdaMetals(high_T, high_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", high_T, high_rho, lambda_6);
  //  double lambda_7 = LambdaMetals(high_T, low_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", high_T, low_rho, lambda_7);
  //  double lambda_8 = LambdaMetals(low_T, high_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", low_T, high_rho, lambda_8);
  //  double lambda_9 = LambdaMetals(T_0, rho_0, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", T_0, rho_0, lambda_9);
  //  double lambda_10 = LambdaMetals(neg_T, neg_rho, MetalT.lambda_m, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("Interpolated cooling for T = %0.3e, nh = %0.3e: Lambda_z = %0.3e\n", neg_T, neg_rho, lambda_10);

  //  double ne1 = neSolar(min_T, max_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", min_T, max_rho, ne1);
  //  double ne2 = neSolar(max_T, min_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", max_T, min_rho, ne2);
  //  double ne3 = neSolar(min_T, min_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", min_T, min_rho, ne3);
  //  double ne4 = neSolar(max_T, max_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", max_T, max_rho, ne4);
  //  double ne5 = neSolar(low_T, low_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", low_T, low_rho, ne5);
  //  double ne6 = neSolar(high_T, high_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", high_T, high_rho, ne6);  
  //  double ne7 = neSolar(high_T, low_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", high_T, low_rho, ne7);  
  //  double ne8 = neSolar(low_T, high_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", low_T, high_rho, ne8);  
  //  double ne9 = neSolar(T_0, rho_0, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", T_0, rho_0, ne9);  
  //  double ne10 = neSolar(neg_T, neg_rho, MetalT.ne_solar, MetalT.temperature_bins, MetalT.nh_bins, MetalT.rows_cooling, MetalT.columns_cooling);
  //  printf("ne solar for T = %0.3e, nh = %0.3e: ne = %0.3e\n", neg_T, neg_rho, ne10);  

//   struct metallicity metallicity;

//   metallicity_init(&SphP[0].Metallicity[i], &SphP[0].pMetallicity[i], SCALAR_TYPE_PASSIVE);

//   /* save type and relative address */
//   metallicity.type        = type;
//   metallicity.offset      = ((char *)addr) - ((char *)&SphP[0]);
//   metallicity.offset_mass = ((char *)addr_mass) - ((char *)&SphP[0]);

// #else  /* #ifdef MAXSCALARS */
//   return -1;
// #endif /* #ifdef MAXSCALARS #else */

   #endif /* #ifdef METALLIC_COOLING */
   mpi_printf("GFM_COOLING: time, time begin = %le\t%le\n", All.Time, All.TimeBegin);
   All.Time = All.TimeBegin;
   set_cosmo_factors_for_current_time();
 
   IonizeParams();
 
 


 }
 
 /*! \brief Apply the isochoric cooling to all the active gas cells.
  *
  *  \return void
  */
 void cooling_only(void) /* normal cooling routine when star formation is disabled */
 {
   int idx, i;
 
   CPU_Step[CPU_MISC] += measure_time();

   for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
     {
       i = TimeBinsHydro.ActiveParticleList[idx];
       if(i >= 0)
         {
           if(P[i].Mass == 0 && P[i].ID == 0)
             continue; /* skip cells that have been swallowed or eliminated */

#ifdef DEBUG_FIND_CELL
  if (P[i].ID == 16713972) // whatever the crash id is 
  {
    printf("Particle ID = %d found in [Node %d] in [cooling.c] at [cooling_only()] before calling [cool_cell()]: SphP[%d].Utherm=%g, SphP[%d].Density=%g, P[%d].Pos = [%f, %f, %f], time = %f\n", 
      16713972, ThisTask, i, SphP[i].Utherm, i, SphP[i].Density, i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], All.Time);
  }
#endif

/* If enabled, we are restricting cooling to the innermost grid. */
#if defined(NOCOOL_BACKGROUND_GRID) 
           double rad_i = P[i].Pos[0] - 0.5*All.BoxSize - All.GlobalDisplacementVector[0]; 
           double rad_j = P[i].Pos[1] - 0.5*All.BoxSize - All.GlobalDisplacementVector[1];
           double rad_k = P[i].Pos[2] - 0.5*All.BoxSize - All.GlobalDisplacementVector[2];
           if (fabs(rad_i) <= All.BoxLimit && fabs(rad_j ) <= All.BoxLimit && fabs(rad_k) <= All.BoxLimit)
             cool_cell(i);
#else 
             cool_cell(i); // cool cell regardless of location
#endif /* NOCOOL_BACKGROUND_GRID*/
         }
     }
   CPU_Step[CPU_COOLINGSFR] += measure_time();
 }
 
 /*! \brief Apply the isochoric cooling to a given gas cell.
  *
  *  This function applies the normal isochoric cooling to a single gas cell.
  *  Once the cooling has been applied according to one of the cooling models
  *  implemented, the internal energy per unit mass, the total energy and the
  *  pressure of the cell are updated.
  *
  *  \param[in] i Index of the gas cell to which cooling is applied.
  *
  *  \return void
  */
 void cool_cell(int i)
 {
#ifdef DEBUG_FIND_CELL
   if (P[i].ID == 16713972) // whatever the crash id is 
   {
    printf("Particle ID = %d found in [Node %d] in [cooling.c] at [cool_cell()] before applying isochoric cooling: SphP[%d].Utherm=%g, SphP[%d].Density=%g, P[%d].Pos = [%f, %f, %f], time = %f\n", 
      16713972, ThisTask, i, SphP[i].Utherm, i, SphP[i].Density, i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], All.Time);
   }
#endif
   double dt, dtime, ne = 1;
   double unew, dens, dtcool;
   double met;
   dens = SphP[i].Density;
 
   dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
 
   dtime = All.cf_atime * dt / All.cf_time_hubble_a;

#ifdef METALLIC_COOLING
   if (SphP[i].PScalars[0] < All.InitBackgroundMetallicity/2)
      SphP[i].PScalars[0] = All.InitBackgroundMetallicity/2;
    // SphP[i].PScalars[0] = 0; //  In case the metallicity is less than 0, then we set it to 0
   met = SphP[i].PScalars[0];
   // Change to terminate instead
  //  met = SphP[i].Metallicity;
#else 
   met = 0;
#endif /* #ifdef METALLIC_COOLING */

   dtcool = dtime;

  //  if (SphP[i].Utherm <= All.MinEgySpec/1e6)
  //  {
  //   printf("All.MinEgySpec=%f\n", All.MinEgySpec);
  //   double rad_i = P[i].Pos[0] - 0.5*All.BoxSize - All.GlobalDisplacementVector[0];
  //   double rad_j = P[i].Pos[1] - 0.5*All.BoxSize - All.GlobalDisplacementVector[1];
  //   double rad_k = P[i].Pos[2] - 0.5*All.BoxSize - All.GlobalDisplacementVector[2];
  
  //   terminate("SphP[%d].Utherm = 0,  P[%d].ID = %d, mass = %f, density = %f, number density(in cm^{-3}) = %f, Radial Coordinate=%f|%f|%f, time=%f", i, i, P[i].ID, P[i].Mass, SphP[i].Density, SphP[i].Density*All.UnitDensity_in_cgs/PROTONMASS, rad_i, rad_j, rad_k, All.Time);     
  //  }

   ne         = SphP[i].Ne; /* electron abundance (gives ionization state and mean molecular weight) */
   unew       = DoCooling(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne, met);
   SphP[i].Ne = ne;
 
   if(unew < 0)
     terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);
 
   double du = unew - SphP[i].Utherm;
 
   if(unew < All.MinEgySpec)
     du = All.MinEgySpec - SphP[i].Utherm;
 
   SphP[i].Utherm += du;
   SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
 
 #ifdef OUTPUT_COOLHEAT
   if(dtime > 0)
     SphP[i].CoolHeat = du * P[i].Mass / dtime;
 #endif /* #ifdef OUTPUT_COOLHEAT */
 
   set_pressure_of_cell(i);
 }
 
 #ifdef METALLIC_COOLING

   /*! \brief Read and extracts values from a PIE metallic cooling HDF5 file.
   *
   *  \param[in] file_name Name of file that contains the tabulated parameters.
   *
   *  \return void */
 void ReadMetallicParams(const char *file_name)
 {
   hsize_t *Tdims, *rhodims;
   hsize_t *gdims_metal;
   hsize_t *nedims;
 
   // getting datasets from the hdf5 file
   hid_t file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT); // open up the file 
   hid_t temp_bin_id = H5Dopen(file_id, "/Total_Metals/Temperature_bins"); // get the temperature dataset
   hid_t nh_bin_id = H5Dopen(file_id, "/Total_Metals/Hydrogen_density_bins"); // get the hydrogen density dataset
   hid_t gm_id_metal = H5Dopen(file_id, "/Total_Metals/Net_cooling"); // get the cooling function data set
   hid_t ne_id = H5Dopen(file_id, "/Solar/Electron_density_over_n_h"); // get the solar electron abundnace 
  //  printf("hid_t id %e\n", ne_id);

   hid_t tbin_space = H5Dget_space(temp_bin_id); // get the space id
   int T_rank = H5Sget_simple_extent_ndims(tbin_space); // get the rank -> 1 for rho or T, 2 for gamma
   Tdims = malloc(T_rank * sizeof(hsize_t)); // dynamicalloy allocate size for T 
   H5Sget_simple_extent_dims(tbin_space, Tdims, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT.temperature_bins = malloc(Tdims[0] * sizeof(double)); // allocate space for the temperature bins
   H5Dread(temp_bin_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalT.temperature_bins); // read and get the temperature bins
   // double temp_min = MetalT.temperature_bins[0]; 
   // double temp_average = 0;
   // double temp_max = MetalT.temperature_bins[0];
   // for (size_t i = 0; i < Tdims[0]; i++)
   // {
   //   if (MetalT.temperature_bins[i] < temp_min)
   //     temp_min = MetalT.temperature_bins[i]; 
   //   if (MetalT.temperature_bins[i] > temp_max)
   //     temp_max = MetalT.temperature_bins[i]; 
   //   temp_average += MetalT.temperature_bins[i];
   // }
   // printf("final minimum temperature is: %e\n", temp_min);
   // printf("final maximum temperature is: %e\n", temp_max);
   // printf("Average temperature bin value is: %e\n", temp_average/Tdims[0]);
 
   hid_t nhbin_space = H5Dget_space(nh_bin_id);  // get the space id
   int rho_rank = H5Sget_simple_extent_ndims(nhbin_space); // get the rank -> 1 for rho or T, 2 for gamma
   rhodims = malloc(rho_rank * sizeof(hsize_t)); // dynamically allocate size for rho
   H5Sget_simple_extent_dims(nhbin_space, rhodims, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT.nh_bins = malloc(rhodims[0] * sizeof(double)); // allocate space for the rho bins
   H5Dread(nh_bin_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalT.nh_bins);  // read and get the hydrogen density bins
   // double rho_min = MetalT.nh_bins[0]; 
   // double rho_average = 0;
   // double rho_max = MetalT.nh_bins[0];
   // for (size_t i = 0; i < rhodims[0]; i++)
   // {
   //   if (MetalT.nh_bins[i] < rho_min)
   //     rho_min = MetalT.nh_bins[i]; 
   //   if (MetalT.nh_bins[i] > rho_max)
   //     rho_max = MetalT.nh_bins[i]; 
 
   //   rho_average += MetalT.nh_bins[i];
   // }
   // printf("final minimum rho is: %e\n", rho_min);
   // printf("final maximum rho is: %e\n", rho_max);
   // printf("Average rho bin value is: %e\n", rho_average/rhodims[0]);
 
   hid_t gm_space = H5Dget_space(gm_id_metal); // get the space id
   int gm_rank = H5Sget_simple_extent_ndims(gm_space); // get the rank -> 1 for rho or T, 2 for gamma
   gdims_metal = malloc(gm_rank * sizeof(hsize_t)); // dynamically allocate size for gamma
   H5Sget_simple_extent_dims(gm_space, gdims_metal, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT.lambda_m = malloc(gdims_metal[0] * gdims_metal[1] * sizeof(double)); // allocate space for the gamma bins
   H5Dread(gm_id_metal, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalT.lambda_m); // read and get the gamma values. Note that this is actually a 1d array with row major format
 
   MetalT.rows_cooling = gdims_metal[0];
   MetalT.columns_cooling = gdims_metal[1];

   hid_t ne_space = H5Dget_space(ne_id); // get the space id
   int ne_rank = H5Sget_simple_extent_ndims(ne_space); // get the rank -> 1 for rho or T, 2 for gamma or ne 
   nedims = malloc(ne_rank * sizeof(hsize_t)); // dynamically allocate size for gamma
   H5Sget_simple_extent_dims(ne_space, nedims, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT.ne_solar = malloc(nedims[0] * nedims[1] * sizeof(double)); // allocate space for the gamma bins
   H5Dread(ne_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalT.ne_solar); // read and get the gamma values. Note that this is actually a 1d array with row major format

   free(Tdims); 
   free(rhodims);
   free(gdims_metal); 
   free(nedims);

   H5Sclose(tbin_space);
   H5Sclose(nhbin_space);
   H5Sclose(gm_space);
   H5Sclose(ne_space);

   H5Dclose(temp_bin_id);
   H5Dclose(nh_bin_id);
   H5Dclose(gm_id_metal);
   H5Dclose(ne_id);

   H5Fclose(file_id);
 }
 
  /*! \brief Returns the interpolated metallic cooling function of a gas cell using a PIE cooling table.
  *
  *  This function takes values from a sorted table and uses a bilinear interpolate to get the metallic cooling rate for a cell. 
  *  It is assumed that the table is small and any performance differences between search algorithms is negligible.
  * 
  *  \param[in] T The temperature of the cell
  *  \param[in] nh The number density of the cell
  *  \param[in] lambda The cooling rates for a given range of temperatures and densities, formated as a 1d array of length (T_len + nh_len)
  *  \param[in] T_bins The temperature bins
  *  \param[in] nh_bins The hydrogen number density bins 
  *  \param[in] T_len The number of the temperature bins
  *  \return interpolated metallic cooling rate given cell temperature and density.
  */
 double LambdaMetals(double T, double nh, double lambda[], double T_bins[], double nh_bins[], size_t T_len, size_t nh_len)
 {

   // 1. Determine the temperature and nh bins and points
   double Tp, nhp; // the Tp points and the rhop points of our simulation
   size_t Th_index = 0, nhh_index = 0;
   // 1.1.a Check if T is outside T_bins and set Tp and Th_index accordingly
   if (T >= T_bins[T_len - 1]) // if T is higher than the highest value, use the last two bins
   {
    Th_index = T_len - 1;
    Tp = T_bins[T_len - 1]; // use the maximum temperature in the table
   }
   else if (T <= T_bins[0])  // if T is lower than the lowest value, use the first two bins
   {
    Th_index = 1;
    Tp = T_bins[0]; // Use the lowest temperature in the table
   }
   // 1.1.b Apply a linear search to search through our temperature bins 
   else 
   {
    for (size_t i = 0; i < T_len; i++)
    {
      Th_index = i;
      if (T_bins[i] > T)
        break;
    }
    Tp = T; // use the cell temperature
   }

   // 1.2.a Check if nh is outside nh_bins and set nhp and nhh_index accordingly
   if (nh >= nh_bins[nh_len - 1]) // if rho is higher than the highest value, use the last two bins
   {
    nhh_index = nh_len - 1;
    nhp = nh_bins[nh_len - 1]; // use the maximum density in the table
   }
   else if (nh <= nh_bins[0]) // if rho is lower than the lowest value, use the first two binss
   {
    nhh_index = 1;
    nhp = nh_bins[0]; // use the lowest density in the table
   }
   // 1.2.b Apply a linear search to our number density bins
   else 
   {
    for (size_t i = 0; i < nh_len; i++)
    {
      nhh_index = i;
      if (nh_bins[i] > nh)
        break;
    }
    nhp = nh; // use the cell number density
   }

   // 2. Get the values for specific bins
   size_t Tl_index = Th_index - 1;
   double T_b1 = T_bins[Tl_index], T_b2 = T_bins[Th_index]; 
   size_t nhl_index = nhh_index - 1;
   double nh_b1 = nh_bins[nhl_index], nh_b2 = nh_bins[nhh_index]; 

   // 3. Perform a bilinear interpolation to get the cooling rate for the disk
   double g11 = lambda[Tl_index*nh_len + nhl_index], g12 = lambda[Tl_index*nh_len + nhh_index], g21 = lambda[Th_index*nh_len + nhl_index], g22 = lambda[Th_index*nh_len + nhh_index];
   double gt1 = (T_b2 - Tp)/(T_b2 - T_b1)*g11 + (Tp - T_b1)/(T_b2 - T_b1)*g21;
   double gt2 =  (T_b2 - Tp)/(T_b2 - T_b1)*g12 + (Tp - T_b1)/(T_b2 - T_b1)*g22;

   return ((nh_b2 - nhp)/(nh_b2 - nh_b1)*gt1 + (nhp - nh_b1)/(nh_b2 - nh_b1)*gt2); 
 }

  /*! \brief Returns the interpolated solar electron abundances for a cell using a PIE cooling table.
  *
  *  This function takes read values from a sorted table and uses a nearest neighbors approach to interpolate the electron abundances for a given cell.
  *  It is assumed that the table is small and any performance differences between search algorithms is negligible.
  *
  *  \param[in] T The temperature of the cell
  *  \param[in] nh The number density of the cell
  *  \param[in] ne_solar The electron abundances for a given range of temperatures and densities, formated as a 1d array of length (T_len + nh_len)
  *  \param[in] T_bins The temperature bins
  *  \param[in] nh_bins The hydrogen number density bins 
  *  \param[in] T_len The number of temperature bins
  *  \param[in] nh_len The number of hydrogen density bins
  *  \return solar electron abundance given cell temperature and density.
  */
 double neSolar(double T, double nh, double ne_solar[], double T_bins[], double nh_bins[], size_t T_len, size_t nh_len)
 {
   double dT = abs(T - T_bins[0]);
   double dnh = abs(nh - nh_bins[0]);
   size_t Tn_index = 0, nhn_index = 0;
   for (size_t i = 1; i < T_len; i++)
   {
    if (abs(T - T_bins[i]) < dT)
    {
      dT = abs(T - T_bins[i]);
      Tn_index = i;
    }
   }
   for (size_t i = 1; i < nh_len; i++)
   {
    if (abs(nh - nh_bins[i]) < dnh)
    {
      dnh = abs(nh - nh_bins[i]);
      nhn_index = i;
    }
   }

   return ne_solar[Tn_index*nh_len + nhn_index]; 
 }


   /*! \brief Read and extracts values from a CIE metallic cooling HDF5 file.
   *
   *  \param[in] file_name Name of file that contains the tabulated parameters.
   *
   *  \return void */
 void ReadMetallicParams_CIE(const char *file_name)
 {
   hsize_t *Tdims;
   hsize_t *gdims_metal;
   hsize_t *nedims;
 
   hid_t file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT); // open up the file 
   hid_t temp_bin_id = H5Dopen(file_id, "/Total_Metals/Temperature_bins"); // get the temperature dataset
   hid_t gm_id_metal = H5Dopen(file_id, "/Total_Metals/Net_cooling"); // get the cooling function data set
   hid_t ne_id = H5Dopen(file_id, "/Solar/Electron_density_over_n_h"); // get the solar electron abundnace 
   // /Total_Metals/Net_cooling
   hid_t tbin_space = H5Dget_space(temp_bin_id); // get the space id
   int T_rank = H5Sget_simple_extent_ndims(tbin_space); // get the rank -> 1 for rho or T, 2 for gamma
   Tdims = malloc(T_rank * sizeof(hsize_t)); // dynamicalloy allocate size for T 
   printf("Tdims: %f\n", Tdims);
   H5Sget_simple_extent_dims(tbin_space, Tdims, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT_CIE.temperature_bins = malloc(Tdims[0] * sizeof(double)); // allocate space for the temperature bins
   H5Dread(temp_bin_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MetalT_CIE.temperature_bins); // read and get the temperature bins
 
   hid_t gm_space = H5Dget_space(gm_id_metal); // get the space id
   int gm_rank = H5Sget_simple_extent_ndims(gm_space); // get the rank -> 1 for rho or T, 2 for gamma
   gdims_metal = malloc(gm_rank * sizeof(hsize_t)); // dynamically allocate size for gamma
   H5Sget_simple_extent_dims(gm_space, gdims_metal, NULL); // Retrieves dataspace dimension size and maximum size
   MetalT_CIE.lambda_m = malloc(gdims_metal[0]* sizeof(double)); // allocate space for the gamma bins
   printf("gdims_metal: %f\n", gdims_metal);

   H5Dread(gm_id_metal, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,  MetalT_CIE.lambda_m); // read and get the gamma values. Note that this is actually a 1d array with row major format
 
    MetalT_CIE.rows_cooling = gdims_metal[0];
  //   MetalT_CIE.columns_cooling = gdims_metal[1];

   hid_t ne_space = H5Dget_space(ne_id); // get the space id
   int ne_rank = H5Sget_simple_extent_ndims(ne_space); // get the rank -> 1 for rho or T, 2 for gamma or ne 
   nedims = malloc(ne_rank * sizeof(hsize_t)); // dynamically allocate size for gamma
   H5Sget_simple_extent_dims(ne_space, nedims, NULL); // Retrieves dataspace dimension size and maximum size
    MetalT_CIE.ne_solar = malloc(nedims[0] * sizeof(double)); // allocate space for the gamma bins
   H5Dread(ne_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,  MetalT_CIE.ne_solar); // read and get the gamma values. Note that this is actually a 1d array with row major format

   free(Tdims); 
   free(gdims_metal); 
   free(nedims);

   H5Sclose(tbin_space);
   H5Sclose(gm_space);
   H5Sclose(ne_space);

   H5Dclose(temp_bin_id);
   H5Dclose(gm_id_metal);
   H5Dclose(ne_id);

   H5Fclose(file_id);
 }
 
  /*! \brief Returns the interpolated metallic cooling function of a gas cell using a CIE cooling table.
  *
  *  This function takes read values from a sorted table and uses a linear interpolate to get the metallic cooling rate for a cell. 
  *  It is assumed that the table is small and any performance differences between search algorithms is negligible.
  * 
  *  \param[in] T The temperature of the cell
  *  \param[in] lambda The cooling rates for a given range of temperatures and densities, formated as a 1d array of length (T_len + nh_len)
  *  \param[in] T_bins The temperature bins
  *  \param[in] T_len The number of the temperature bins
  *  \return interpolated cooling rate given cell temperature.
  */
 double LambdaMetals_CIE(double T, double lambda[], double T_bins[], size_t T_len)
 {
   // 1. Determine the temperature  and points
   double Tp;
   size_t Th_index = 0;
   // 1.1.a Check if T is outside T_bins and set Tp and Th_index accordingly
   if (T >= T_bins[T_len - 1]) // if T is higher than the highest value, use the last two bins
   {
    Th_index = T_len - 1;
    Tp = T_bins[T_len - 1]; // use the maximum temperature in the table
   }
   else if (T <= T_bins[0])  // if T is lower than the lowest value, use the first two bins
   {
    Th_index = 1;
    Tp = T_bins[0]; // Use the lowest temperature in the table
   }
   // 1.1.b Apply a linear search to search through our temperature bins 
   else 
   {
    for (size_t i = 0; i < T_len; i++)
    {
      Th_index = i;
      if (T_bins[i] > T)
        break;
    }
    Tp = T; // use the cell temperature
   }

   // 2. Get the values for specific bins
   size_t Tl_index = Th_index - 1;
   double T1 = T_bins[Tl_index], T2 = T_bins[Th_index]; 

   // 3. Perform a linear interpolation to get the cooling rate for the disk
   // f(x) = f(x1)*(x2 -x)/(x2 - x1) + f(y2)(x-x1)/(x2-x1)
   double g1 = lambda[Tl_index];
   double g2 = lambda[Th_index];
   return ((T2 - Tp)/(T2 - T1))*g1 + ((Tp - T1)/(T2 - T1))*g2;
 }

  /*! \brief Returns the interpolated solar electron abundances for a cell using a CIE cooling table.
  *
  *  This function takes read values from a sorted table and uses a nearest neighbors approach to interpolate the electron abundances for a given cell.
  *  It is assumed that the table is small and any performance differences between search algorithms is negligible.
  *
  *  \param[in] T The temperature of the cell
  *  \param[in] ne_solar The electron abundances for a given range of temperatures and densities, formated as a 1d array of length (T_len + nh_len)
  *  \param[in] T_bins The temperature bins
  *  \param[in] T_len The number of the temperature bins
  *  \return solar electron abundance given cell temperature
  */
 double neSolar_CIE(double T, double ne_solar[], double T_bins[], size_t T_len)
 {
   double dT = abs(T - T_bins[0]);
   size_t Tn_index = 0;
   for (size_t i = 1; i < T_len; i++)
   {
    if (abs(T - T_bins[i]) < dT)
    {
      dT = abs(T - T_bins[i]);
      Tn_index = i;
    }
   }
   return ne_solar[Tn_index]; 
 }
 #endif /* #ifdef METALLIC_COOLING */

#endif /* #ifdef COOLING */
 