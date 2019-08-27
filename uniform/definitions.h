#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 3
#define  USER_DEF_PARAMETERS     15
#define  USER_DEF_CONSTANTS      9

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          YES
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA_GAS               0
#define  CLOUD_RADIUS            1
#define  CLOUD_CENTREX1          2
#define  CLOUD_CENTREX2          3
#define  CLOUD_CENTREX3          4
#define  RHO_W                   5
#define  RHO_C                   6
#define  PRS_C                   7
#define  VX1_W                   8
#define  VX1_C                   9
#define  BETA_W                  10
#define  BETA_C                  11
#define  DENSITY_STEEPNESS       12
#define  CORE_RADIUS             13
#define  CUT_RADIUS              14

/* -- user-defined symbolic constants -- */

#define  UNIT_LENGTH             1000*CONST_pc
#define  FIELD_NONE              0
#define  FIELD_ORDERED           1
#define  FIELD_TANGLED           2
#define  FIELD_TRANSVERSE        1
#define  FIELD_PARALLEL          2
#define  FIELD_OBLIQUE           3
#define  CLOUD_FIELD             FIELD_TRANSVERSE
#define  BG_FIELD                FIELD_TRANSVERSE

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          MULTID
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   DEFAULT
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
