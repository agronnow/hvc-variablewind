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
#define ISOTHERMAL_WIND

#ifdef VARIABLE_INFLOW
double g_tabulardistance[1024];
double g_tabulardensity[1024];
double g_tabularbfield[1024];
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
	  if (i > 1023) {return -2;}
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
  #endif

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
#ifdef OUTP_BOV
  double dt = 0.1;
  int filenumber = g_time/dt;
  if (prank == 0)
  {
    char* dir;
    char filename[512];
    FILE* BOVout;
    dir = GetOutputDir();
    Output vnames;
    Output* varnames = &vnames;
    varnames->var_name = ARRAY_2D(64,128,char);
    SetDefaultVarNames(varnames);

    int nv = NVAR;
    #ifdef STAGGERED_MHD
     D_EXPAND(
       varnames->var_name[nv]   = "bx1s";
       nv++; ,

       varnames->var_name[nv]   = "bx2s"; 
       nv++; ,

       varnames->var_name[nv]   = "bx3s"; 
       nv++;
     )
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     #if DIMENSIONS == 3   
      varnames->var_name[nv]   = "Ax1";
      nv++;
      varnames->var_name[nv]   = "Ax2";
      nv++;
     #endif
     varnames->var_name[nv]   = "Ax3";
     nv++;
    #endif

    int i;
    for (i = 0; i < nv; i++)
    {
      sprintf(filename, "%s/%s.%04d%s", dir, varnames->var_name[i], filenumber, ".bov");
      BOVout = fopen(filename, "w");
      fprintf(BOVout, "TIME: %f\nDATA_FILE: %s.%04d%s\nDATA_SIZE: %d %d %d\nDATA_FORMAT: DOUBLE\n"
	      "VARIABLE: %s\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n"
	      "BRICK_ORIGIN: %f %f %f\nBRICK_SIZE: %f %f %f",
	      g_time, varnames->var_name[i], filenumber, ".dbl", grid[IDIR].np_int_glob, grid[JDIR].np_int_glob, grid[KDIR].np_int_glob,
	      varnames->var_name[i], g_inputParam[CLOUD_CENTREX1], g_inputParam[CLOUD_CENTREX2], g_inputParam[CLOUD_CENTREX3],
	      g_domEnd[IDIR]-g_domBeg[IDIR], g_domEnd[JDIR]-g_domBeg[JDIR], g_domEnd[KDIR]-g_domBeg[KDIR]);
      fclose(BOVout);
    }

//    filenumber++;
  }
#endif //OUTP_BOV


#ifdef OUTP_CENTER_OF_MASS
  int i, j, k;
    int var = NVAR;
    //  print1("var %f", d->Vc[TRC][middlez][middley][middlex+1]);

    double LocalCenterOfMassX = 0.0;
    double LocalMass = 0.0;
    double GlobalCenterOfMassX = 0.0;
    double GlobalMass = 0.0;
    //double dV = grid[IDIR].dx[0] * grid[JDIR].dx[0] * grid[KDIR].dx[0];

    //    print1("dV %f", grid[IDIR].dx[0] * UNIT_LENGTH);
    //print("rank %i, gbeg %i, gend %i ,ghost %i, globbeg %f globend %f ",prank, grid[IDIR].beg,grid[IDIR].end,grid[IDIR].nghost,grid[IDIR].x_glob[grid[IDIR].beg],grid[IDIR].x_glob[grid[IDIR].end]);
    //    print("rank %i, x1 %f, x2 %f ",prank,(grid[IDIR].x_glob[(grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[0]/2.0), (grid[IDIR].x_glob[IEND + (grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[IEND]/2.0));

    TOT_LOOP(k,j,i)
    {
      LocalMass += d->Vc[RHO][k][j][i];// * dV; //For a uniform grid the volume term is constant cancels out when computing the global center of mass
      int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
      LocalCenterOfMassX += (grid[IDIR].x_glob[GlobalIndex] + grid[IDIR].dx[i]/2.0) * d->Vc[RHO][k][j][i];// * dV;
    }

    /* TODO: MPI STUFF SO THAT THIS WORKS IN PARALLEL
             Hopefully this won't slow the processes down too much. If it does, maybe find a way to only do this ocassionaly rather than at every timestep*/\

    MPI_Allreduce(&LocalCenterOfMassX, &GlobalCenterOfMassX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalMass, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    GlobalCenterOfMassX /= GlobalMass;

    double maxrho = -1.0;
    double maxrhoX = 0.0;
    TOT_LOOP(k, j, i)
    {
      if (d->Vc[RHO][k][j][i] > maxrho)
      {
	int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
	maxrho = d->Vc[RHO][k][j][i];
	maxrhoX = grid[IDIR].x_glob[GlobalIndex];
      }
    }

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    double *GlobalmaxrhoX = (double *)malloc(sizeof(double) * numprocs);
    double *Globalmaxrho = (double *)malloc(sizeof(double) * numprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&maxrhoX, 1, MPI_DOUBLE, GlobalmaxrhoX, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&maxrho, 1, MPI_DOUBLE, Globalmaxrho, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    int maxelem = -1.0;
    int elem;
    for (elem=0; elem < numprocs; elem++)
    {
      print("Glob %f ", Globalmaxrho[elem]);
	 if (Globalmaxrho[elem] > Globalmaxrho[maxelem])
	 {
	    maxelem = elem;
	 }
    }

    maxrhoX = GlobalmaxrhoX[maxelem];
    
    free(GlobalmaxrhoX);
    free(Globalmaxrho);

    if (prank == 0)
      {
	char *dir, fname[512];
	FILE *fp;
	dir = GetOutputDir();
	sprintf(fname, "%s/centerofmassx.txt", dir);
	fp = fopen(fname, "a");
	fprintf(fp, "%12.6e  %12.6e %12.6e\n", g_time, GlobalCenterOfMassX, maxrhoX);
	fclose(fp);
      }
#endif //OUTP_CENTER_OF_MASS

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
  double RhoCloud = g_inputParam[RHO_C];
#ifdef PRS_FROM_TEMP
  double PrsWind = RhoCloud*TCloud/(KELVIN*muCloud);
#else
  double PrsWind = g_inputParam[PRS_C];
#endif

#ifdef VARIABLE_INFLOW
  double PrsCloud = PrsWind;
  static double InitRhoWind = 0;
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
        if (dist > init_dist - g_tabulardistance[0]) {g_curidx = 0;}
	else
	{
  	  while (init_dist - dist < g_tabulardistance[g_curidx])
	  {
	    g_curidx--;
	  }
	}
	if (InitRhoWind == 0) {InitRhoWind = g_tabulardensity[g_curidx];}
        if (niter == 100000) {print1("Info: %f %f %f %f\n", 50.0-dist, g_tabulardistance[g_curidx], g_tabulardensity[g_curidx], g_tabularbfield[g_curidx]);niter=0;}
        niter++;
	d->Vc[RHO][k][j][i] = g_tabulardensity[g_curidx];
	#endif
	double B0 = 0;
        #if PHYSICS == MHD
	#ifdef VARIABLE_INFLOW
	B0 = g_tabularbfield[g_curidx] * MICROGAUSS;
	#else
	B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
	#endif
        #endif
	d->Vc[VX1][k][j][i] = VelWind;
        #ifdef ISOTHERMAL_WIND
	d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*PrsCloud/InitRhoWind;
	#endif
	#if PHYSICS == MHD
	EXPAND(d->Vc[BX1][k][j][i] = 0; ,
	       d->Vc[BX2][k][j][i] = B0; ,
	       d->Vc[BX3][k][j][i] = 0;)
//	print("B0: %f\n",B0);
	#endif
      }
    }
    else if (box->vpos == X1FACE) {}
    else if (box->vpos == X2FACE)
    {
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
	#ifdef VARIABLE_INFLOW
	d->Vs[BX2s][k][j][i] = g_tabularbfield[g_curidx] * MICROGAUSS;
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
	d->Vs[BX3s][k][j][i] = 0;//g_tabularbfield[g_curidx]/sqrt3 * MICROGAUSS;
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
