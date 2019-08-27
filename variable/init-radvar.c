/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author Asger Gronnow (asger.gronnow@usyd.edu.au)
  \date   Nov 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#define PRS_FROM_TEMP
//#define VELOCITY_MACH
#define TRC_CLOUD TRC
#define TRC_CORE TRC+1
#define TRC_ENV TRC+2

#define MICROGAUSS 2.181
#define VARIABLE_INFLOW

#ifdef VARIABLE_INFLOW
double g_tabulardistance[2048];
double g_tabulardensity[2048];
double g_tabularbfield[2048];
int g_curidx = 0;
#endif

/* ********************************************************************* */
int ReadTable(char TableFile[256])
/*!
 * Read table of z-distance, halo density and magnetic field strength
 * used when calculating variable inflow.
 *
 * \param [in] TableFile     Path to file containing table with distance
 *                           in kpc in the first column, halo number
 *                           density per cubic cm in the second column
 *                           and halo magnetic field magnitude in
 *                           microGauss in the third column.
 *
 * \return 0 on success or one of the following negative integers in
 *         case of an error:
 *         -1 Table file could not be opened
 *         -2 Table has too many data points (> 1024)
 *********************************************************************** */
{
  int i = 0;
  FILE* table = fopen(TableFile, "r");
  if (table == NULL) {return -1;}
  else
    {
      char line[256];
      double dist;
      double dens;
      double bfield;
      while (fgets(line, sizeof(line), table) != NULL) 
	{
	  if (line[0] != '#')
	    {
	      sscanf(line, "%lf %lf %lf", &dist, &dens, &bfield);
	      g_tabulardistance[i] = dist;
	      g_tabulardensity[i] = dens;
	      g_tabularbfield[i] = bfield;
	      i++;
	    }
	  if (i > 2047) {return -2;}
	}
    }
  fclose(table);
  g_curidx = i - 1;
  return 0;
}


double CenterOfMassX(const Data *d, Grid *grid)
{
  int i, j, k;
  int var = NVAR;

  double LocalCenterOfMassX = 0.0;
  double LocalMass = 0.0;
  double GlobalCenterOfMassX = 0.0;
  double GlobalMass = 0.0;

  DOM_LOOP(k,j,i)
  {
    LocalMass += d->Vc[RHO][k][j][i];// * dV; //For a uniform grid the volume term is constant cancels out when computing the global center of mass
    int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
    LocalCenterOfMassX += (grid[IDIR].x_glob[GlobalIndex] + grid[IDIR].dx[i]/2.0) * d->Vc[RHO][k][j][i];// * dV;
  }

  print("XCM: %d %f\n",prank, LocalCenterOfMassX);

  MPI_Allreduce(&LocalCenterOfMassX, &GlobalCenterOfMassX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&LocalMass, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  GlobalCenterOfMassX /= GlobalMass;
  return GlobalCenterOfMassX;
}


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  static int firstcall = 1;
  double CloudRad = g_inputParam[CLOUD_RADIUS];
  double CoreRad = g_inputParam[CORE_RADIUS];
  double CutRad = g_inputParam[CUT_RADIUS];
  double CurRad;
  double TCloud = 1.e3;
  double muCloud = 1.2;
  double RhoWind = g_inputParam[RHO_W];
  double RhoCloud = g_inputParam[RHO_C];
  g_gamma = g_inputParam[GAMMA_GAS];
  #ifdef PRS_FROM_TEMP
  double PrsCloud = RhoCloud*TCloud/(KELVIN*muCloud);
  #else
  double PrsCloud = g_inputParam[PRS_C];
  #endif

  #if PHYSICS == MHD
  double B0 = sqrt(2*PrsCloud/g_inputParam[BETA_W]);
  #endif

  g_smallDensity   = 1.0000e-07;  // included to avoid negative densities at cloud's front
  g_smallPressure  = 1.0000e-07;  // included to avoid negative pressures at cloud's front


  #ifdef VARIABLE_INFLOW
  #ifdef NO_RADIAL_VAR
  if (firstcall)
  {
    int tableret = ReadTable("distdensbfield-sun10.txt");
    if (tableret != 0) {print("Error reading table on proc %d: %d\n",prank,tableret);}
    firstcall = 0;
  }
  #if PHYSICS == MHD
  B0 = g_tabularbfield[g_curidx] * MICROGAUSS;
  #endif
  RhoWind = g_tabulardensity[g_curidx];
  #else
  if (firstcall)
  {
    int tableret = ReadTable("distdens-radial.txt");
    if (tableret != 0) {print("Error reading table on proc %d: %d\n",prank,tableret);}
    firstcall = 0;
  }
  double init_dist = 50.0;
  //Galactocentric coordinates
  double xg = x3 + 4.0;
  double yg = x2;
  double zg = x1 + init_dist + 2.0;
  double R2 = xg*xg + yg*yg;
  double r2 = R2 + zg*zg;
  int idx = g_curidx;
  if (r2 < g_tabulardistance[0]*g_tabulardistance[0]) {idx = 0;}
  else
  {
    while (r2 < g_tabulardistance[idx]*g_tabulardistance[idx])
    {
      idx--;
    }
  }
  RhoWind = g_tabulardensity[idx];
  //  print1("rho, idx, r2, R2: %f %d %f %f\n",RhoWind,idx,r2,R2);
  #endif

  #endif //VARIABLE_INFLOW

  double x = x1 - g_inputParam[CLOUD_CENTREX1];
  double y = x2 - g_inputParam[CLOUD_CENTREX2];
  double z = x3 - g_inputParam[CLOUD_CENTREX3];

  CurRad = EXPAND(x*x , + y*y , + z*z);
  CurRad = sqrt(CurRad);

  v[TRC_CLOUD] = 0;
  v[TRC_CORE] = 0;
  v[TRC_ENV] = 0;


  if (CurRad < CutRad)
  {
    v[RHO] = RhoWind + (RhoCloud - RhoWind)/(1+pow(CurRad/CoreRad,g_inputParam[DENSITY_STEEPNESS]));
    v[PRS] = PrsCloud;
    if (CurRad < CloudRad)
    {
      EXPAND(v[VX1] = g_inputParam[VX1_C]; ,
	     v[VX2] = 0.0; ,
	     v[VX3] = 0.0;)
      v[TRC_CLOUD] = 1.0;
      if (CurRad < CoreRad) {v[TRC_CORE] = 1.0;}
      else {v[TRC_ENV] = 1.0;}
    }
    else
    {
    #ifdef VELOCITY_MACH
      double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsCloud/RhoWind);
    #else
      double VelWind = g_inputParam[VX1_W];
    #endif
      EXPAND(v[VX1] = VelWind; ,
	     v[VX2] = 0.0; ,
	     v[VX3] = 0.0;)

    }
  }
  else
  {
    #ifdef VELOCITY_MACH
    double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsCloud/RhoWind);
    #else
    double VelWind = g_inputParam[VX1_W];
    #endif
    v[RHO] = RhoWind;
    v[PRS] = PrsCloud;
    EXPAND(v[VX1] = VelWind; ,
           v[VX2] = 0.0; ,
           v[VX3] = 0.0;)
  }


#if PHYSICS == MHD
  #if defined(VARIABLE_INFLOW) && !defined(NO_RADIAL_VAR)
  double zprime = (zg - 1.5)/4.0;
  double R = sqrt(R2);
  //Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
  double Bphi = 2.0/(1+zprime*zprime)*(R/4.0)*exp(-(R-4.0)/4.0);
  if (xg < 1.e-16) {xg = 1.e-16;} //Tiny lower limit on x to avoid division by zero
  //Convert from cylindrical to cartesian simulation coords
  double phi = atan2(yg,xg);
  double Bz = -sin(phi)*Bphi; //simulation coordinate, x-component in galactocentric coordinates
  double By = cos(phi)*Bphi;

  v[BX1] = 0;
  v[BX2] = By * MICROGAUSS;
  v[BX3] = Bz * MICROGAUSS;
  #else
  double Bx = 0;
  double By = 0;
  double Bz = 0;
  double Ax = 0;
  double Ay = 0;
  double Az = 0;

  #if BG_FIELD == FIELD_TRANSVERSE /* Background field is B=(0, B0, 0) */
  By = B0;
  Az = -x*B0;
  #elif BG_FIELD == FIELD_PARALLEL /* Background field is B=(B0, 0, 0) */
  Bx = B0;
  Az = y*B0;
  #elif BG_FIELD == FIELD_OBLIQUE  /* Background field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
  Bx = B0/sqrt(3);
  By = B0/sqrt(3);
  Bz = B0/sqrt(3);
  Ay = B0*x/sqrt(3);
  Az = B0*(y-x)/sqrt(3);
  #endif

  v[BX1] = Bx;
  v[BX2] = By;
  v[BX3] = Bz;

  v[AX1] = Ax;
  v[AX2] = Ay;
  v[AX3] = Az;
  #endif //defined(VARIABLE_INFLOW) && !defined(NO_RADIAL_VAR)

#endif // PHYSICS == MHD
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int i, j, k;
  int nv;

  double TCloud = 1.e3;//2.e4;
//  static double sqrt3 = sqrt(3.0);
  double muCloud = 1.2;
  double RhoWind = g_inputParam[RHO_W];
  #ifdef VARIABLE_INFLOW
  RhoWind = g_tabulardensity[g_curidx];
  #endif
  double RhoCloud = g_inputParam[RHO_C];
#ifdef PRS_FROM_TEMP
  double PrsWind = RhoCloud*TCloud/(KELVIN*muCloud);
#else
  double PrsWind = g_inputParam[PRS_C];
#endif
#ifdef VELOCITY_MACH
  double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsWind/RhoWind);
#else
  double VelWind = g_inputParam[VX1_W];
#endif

  static int niter = 0;

  if (side == 0) {    /* -- check solution inside domain -- */
 //   DOM_LOOP(k,j,i){};
    TOT_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < g_smallDensity) {d->Vc[RHO][k][j][i] = g_smallDensity;}
      if (d->Vc[PRS][k][j][i] < g_smallPressure) {d->Vc[PRS][k][j][i] = g_smallPressure;}
      if (d->Vc[TRC][k][j][i] < 0.0) {d->Vc[TRC][k][j][i] = 0.0;}     
      if (d->Vc[TRC+1][k][j][i] < 0.0) {d->Vc[TRC+1][k][j][i] = 0.0;}
      if (d->Vc[TRC+2][k][j][i] < 0.0) {d->Vc[TRC+2][k][j][i] = 0.0;}
      /*      if (grid[IDIR].x[i] < -1.8)
      {
	d->Vc[BX1][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->Vc[BX2][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->Vc[BX3][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
      }*/
    }
  }
  // set wind (inflow) properties
  else if (side == X1_BEG)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];}
	#ifdef VARIABLE_INFLOW
        double dist = VelWind*g_time;
	double init_dist = 50.0;
	#ifdef NO_RADIAL_VAR
        if (dist > init_dist - g_tabulardistance[0]) {g_curidx = 0;}
	else
	{
  	  while (init_dist - dist < g_tabulardistance[g_curidx])
	  {
	    g_curidx--;
	  }
	}
	d->Vc[RHO][k][j][i] = g_tabulardensity[g_curidx];
        if (niter == 100000) {print1("Info: %f %f %f %f\n", 50.0-dist, g_tabulardistance[g_curidx], g_tabulardensity[g_curidx], g_tabularbfield[g_curidx]);niter=0;}
	#else
	//Galactocentric x,y,z coordinates of the injection plane from simulation coordinates
	int k_global = k + (grid[KDIR].beg - grid[KDIR].nghost);
	int j_global = j + (grid[JDIR].beg - grid[JDIR].nghost);
	double x = grid[KDIR].x_glob[k_global] + 4.0;
	double y = grid[JDIR].x_glob[j_global];
	double z = init_dist - dist;
	double R2 = x*x + y*y;
	double r2 = R2 + z*z;
	int idx = g_curidx;
        if (r2 < g_tabulardistance[0]*g_tabulardistance[0]) {idx = 0;}
	else
	{
  	  while (r2 < g_tabulardistance[idx]*g_tabulardistance[idx])
	  {
	    idx--;
	  }
	}
	d->Vc[RHO][k][j][i] = g_tabulardensity[idx];
        if (niter == 100000) {print1("Info: %f %f %f %f\n", 50.0-dist, g_tabulardistance[idx], g_tabulardensity[idx], g_tabularbfield[idx]);niter=0;}
        #endif //NO_RADIAL_VAR
        niter++;
        #endif //VARIABLE_INFLOW
	double B0 = 0;
        #if PHYSICS == MHD
	#ifdef VARIABLE_INFLOW
        #ifdef NO_RADIAL_VAR
	B0 = g_tabularbfield[g_curidx] * MICROGAUSS;
	EXPAND(d->Vc[BX1][k][j][i] = 0; ,
	       d->Vc[BX2][k][j][i] = B0; ,
	       d->Vc[BX3][k][j][i] = 0;)
        #else
	double zprime = (z - 1.5)/4.0;
	double R = sqrt(R2);
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime)*(R/4.0)*exp(-(R-4.0)/4.0);
	if (x < 1e-16) {x = 1e-16;} //Tiny lower limit on x to avoid division by zero
	//Convert from cylindrical to cartesian simulation coords
	double phi = atan2(y,x);
	double Bz = -sin(phi)*Bphi; //simulation coordinate, x-component in galactocentric coordinates
	double By = cos(phi)*Bphi;
	EXPAND(d->Vc[BX1][k][j][i] = 0; ,
	       d->Vc[BX2][k][j][i] = By * MICROGAUSS; ,
	       d->Vc[BX3][k][j][i] = Bz * MICROGAUSS;)
        #endif //NO_RADIAL_VAR
	#else
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
	EXPAND(d->Vc[BX1][k][j][i] = 0; ,
	       d->Vc[BX2][k][j][i] = B0; ,
	       d->Vc[BX3][k][j][i] = 0;)
        #endif //VARIABLE_INFLOW
        #endif //PHYSICS == MHD
	d->Vc[VX1][k][j][i] = VelWind;
      }
    }
    else if (box->vpos == X1FACE) {}
    else if (box->vpos == X2FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	#ifdef VARIABLE_INFLOW
	#ifdef NO_RADIAL_VAR
	d->Vs[BX2s][k][j][i] = g_tabularbfield[g_curidx] * MICROGAUSS;
	#else
        double dist = VelWind*g_time;
	double init_dist = 50.0;
	int k_global = k + (grid[KDIR].beg - grid[KDIR].nghost);
	int j_global = j + (grid[JDIR].beg - grid[JDIR].nghost);
	double x = grid[KDIR].x_glob[k_global] + 4.0;
	double y = grid[JDIR].x_glob[j_global];
	double z = init_dist - dist;
	double R = sqrt(x*x + y*y);
	double zprime = (z - 1.5)/4.0;
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime)*(R/4.0)*exp(-(R-4.0)/4.0);
	if (x < 1e-16) {x = 1e-16;} //Tiny lower limit on x to avoid division by zero
	//Convert from cylindrical to cartesian simulation coords
	double phi = atan2(y,x);
	double By = cos(phi)*Bphi;
	d->Vs[BX2s][k][j][i] = By * MICROGAUSS;
	#endif
	#else
	d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
	#endif
      }
      #endif
    }
    else if (box->vpos == X3FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	#ifdef VARIABLE_INFLOW
	#ifdef NO_RADIAL_VAR
	d->Vs[BX3s][k][j][i] = 0;//g_tabularbfield[g_curidx]/sqrt3 * MICROGAUSS;
	#else
        double dist = VelWind*g_time;
	double init_dist = 50.0;
	int k_global = k + (grid[KDIR].beg - grid[KDIR].nghost);
	int j_global = j + (grid[JDIR].beg - grid[JDIR].nghost);
	double x = grid[KDIR].x_glob[k_global] + 4.0;
	double y = grid[JDIR].x_glob[j_global];
	double z = init_dist - dist;
	double R = sqrt(x*x + y*y);
	double zprime = (z - 1.5)/4.0;
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime)*(R/4.0)*exp(-(R-4.0)/4.0);
	if (x < 1e-16) {x = 1e-16;} //Tiny lower limit on x to avoid division by zero
	//Convert from cylindrical to cartesian simulation coords
	double phi = atan2(y,x);
	double Bz = -sin(phi)*Bphi; //simulation coordinate, x-component in galactocentric coordinates
	d->Vs[BX3s][k][j][i] = Bz * MICROGAUSS;
	#endif
	#else
	d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IBEG];
	#endif
      }
      #endif
    }
  }
  else if (side == X1_END)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];}
	double B0 = 0;
        #if PHYSICS == MHD
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
        #endif
	d->Vc[VX1][k][j][i] = MAX(d->Vc[VX1][k][j][IEND],0.); // prevent inflow
      }
    }
    else if (box->vpos == X1FACE) {}
    else if (box->vpos == X2FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IEND];
      }
      #endif
    }
    else if (box->vpos == X3FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IEND];
      }
      #endif
    }
  }
  else if (side == X2_BEG)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];}
	double B0 = 0;
        #if PHYSICS == MHD
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
        #endif
	d->Vc[VX2][k][j][i] = MIN(d->Vc[VX2][k][JBEG][i],0.); // prevent inflow
      }
    }
    else if (box->vpos == X1FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX1s][k][j][i] = d->Vs[BX1s][k][JBEG][i];
      }
      #endif
    }
    else if (box->vpos == X2FACE) {}
    else if (box->vpos == X3FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][JBEG][i];
      }
      #endif
    }
  }

  else if (side == X2_END)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][k][JEND][i];}
	double B0 = 0;
        #if PHYSICS == MHD
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
        #endif
	d->Vc[VX2][k][j][i] = MAX(d->Vc[VX2][k][JEND][i],0.); // prevent inflow
      }
    }
    else if (box->vpos == X1FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX1s][k][j][i] = d->Vs[BX1s][k][JEND][i];
      }
      #endif
    }
    else if (box->vpos == X2FACE) {}
    else if (box->vpos == X3FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][JEND][i];
      }
      #endif
    }
  }
  else if (side == X3_BEG)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][KBEG][j][i];}
	double B0 = 0;
        #if PHYSICS == MHD
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
        #endif
	d->Vc[VX3][k][j][i] = MIN(d->Vc[VX3][KBEG][j][i],0.); // prevent inflow
      }
    }
    else if (box->vpos == X1FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX1s][k][j][i] = d->Vs[BX1s][KBEG][j][i];
      }
      #endif
    }
    else if (box->vpos == X2FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX2s][k][j][i] = d->Vs[BX2s][KBEG][j][i];
      }
      #endif
    }
    else if (box->vpos == X3FACE) {}
  }
  else if (side == X3_END)
  {
    if (box->vpos == CENTER)
    {
      BOX_LOOP(box,k,j,i)
      {
	for (nv = 0; nv < NVAR; nv++) {d->Vc[nv][k][j][i] = d->Vc[nv][KEND][j][i];}
	double B0 = 0;
        #if PHYSICS == MHD
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
        #endif
	d->Vc[VX3][k][j][i] = MAX(d->Vc[VX3][KEND][j][i],0.); // prevent inflow
      }
    }
    else if (box->vpos == X1FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX1s][k][j][i] = d->Vs[BX1s][KEND][j][i];
      }
      #endif
    }
    else if (box->vpos == X2FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	d->Vs[BX2s][k][j][i] = d->Vs[BX2s][KEND][j][i];
      }
      #endif
    }
    else if (box->vpos == X3FACE) {}
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
