#include <cplex_consts.h>

using namespace std;

map<std::string, int> cplex_parameters;

void populate_cpx_const()
{

/* --------------------------------------------------------------------------
 * File: cplex.h  
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 1988, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

/* CPX_INFBOUND:  Any bound bigger than this is treated as
   infinity */

//cplex_parameters["CPX_INFBOUND"] = 1.0E+20;


cplex_parameters[" CPX_STR_PARAM_MAX"] =   512;

/* Types of parameters */
 
cplex_parameters["CPX_PARAMTYPE_NONE"] =    0;
cplex_parameters["CPX_PARAMTYPE_INT"] =     1;
cplex_parameters["CPX_PARAMTYPE_DOUBLE"] =  2;
cplex_parameters["CPX_PARAMTYPE_STRING"] =  3;



/* Values returned for 'stat' by solution  */

cplex_parameters["CPX_STAT_OPTIMAL"] = 1;

cplex_parameters["CPX_STAT_OPTIMAL"] =                 1;
cplex_parameters["CPX_STAT_UNBOUNDED"] =               2;
cplex_parameters["CPX_STAT_INFEASIBLE"] =              3;
cplex_parameters["CPX_STAT_INForUNBD"] =               4;
cplex_parameters["CPX_STAT_OPTIMAL_INFEAS"] =          5;
cplex_parameters["CPX_STAT_NUM_BEST"] =                6;
cplex_parameters["CPX_STAT_ABORT_IT_LIM"] =           10;
cplex_parameters["CPX_STAT_ABORT_TIME_LIM"] =         11;
cplex_parameters["CPX_STAT_ABORT_OBJ_LIM"] =          12;
cplex_parameters["CPX_STAT_ABORT_USER"] =             13;

cplex_parameters["CPX_STAT_FEASIBLE_RELAXED_SUM"] =   14;

cplex_parameters["CPX_STAT_OPTIMAL_RELAXED_SUM"] =    15;

cplex_parameters["CPX_STAT_FEASIBLE_RELAXED_INF"] =   16;

cplex_parameters["CPX_STAT_OPTIMAL_RELAXED_INF"] =    17;

cplex_parameters["CPX_STAT_FEASIBLE_RELAXED_QUAD"] =  18;

cplex_parameters["CPX_STAT_OPTIMAL_RELAXED_QUAD"] =   19;

cplex_parameters["CPX_STAT_FEASIBLE"] =                23;


/* Solution type return values from CPXsolninfo */
cplex_parameters["CPX_NO_SOLN"] =       0;
cplex_parameters["CPX_BASIC_SOLN"] =    1;
cplex_parameters["CPX_NONBASIC_SOLN"] = 2;
cplex_parameters["CPX_PRIMAL_SOLN"] =   3;




/* Values of presolve 'stats' for columns and rows */
cplex_parameters["CPX_PRECOL_LOW"] =    -1  /* fixed to original lb */;
cplex_parameters["CPX_PRECOL_UP"] =     -2  /* fixed to original ub */;
cplex_parameters["CPX_PRECOL_FIX"] =    -3  /* fixed to some other value */;
cplex_parameters["CPX_PRECOL_AGG"] =    -4  /* aggregated y = a*x + b */;
cplex_parameters["CPX_PRECOL_OTHER"] =  -5;  /* cannot be expressed by a linear combination;
                                              * of active variables in the presolved model
                                              *  -> crushing will fail if it has to touch
                                              *     such a variable
                                              */
 cplex_parameters["CPX_PREROW_RED"] =    -1  /* redundant row removed in presolved model */;
cplex_parameters["CPX_PREROW_AGG"] =    -2  /* used to aggregate a variable */;
cplex_parameters["CPX_PREROW_OTHER"] =  -3;  /* other, for example merge two inequalities
                                              * into a single equation */


/* Error codes */

cplex_parameters["CPXERR_NO_MEMORY"] =             1001;
cplex_parameters["CPXERR_NO_ENVIRONMENT"] =        1002;
cplex_parameters["CPXERR_BAD_ARGUMENT"] =          1003;
cplex_parameters["CPXERR_NULL_POINTER"] =          1004;
cplex_parameters["CPXERR_CALLBACK"] =              1006;
cplex_parameters["CPXERR_NO_PROBLEM"] =            1009;
cplex_parameters["CPXERR_LIMITS_TOO_BIG"] =        1012;
cplex_parameters["CPXERR_BAD_PARAM_NUM"] =         1013;
cplex_parameters["CPXERR_PARAM_TOO_SMALL"] =       1014;
cplex_parameters["CPXERR_PARAM_TOO_BIG"] =         1015;
cplex_parameters["CPXERR_RESTRICTED_VERSION"] =    1016;
cplex_parameters["CPXERR_NOT_FOR_MIP"] =           1017;
cplex_parameters["CPXERR_NOT_FOR_QP"] =            1018;
cplex_parameters["CPXERR_CHILD_OF_CHILD"] =        1019;
cplex_parameters["CPXERR_TOO_MANY_THREADS"] =      1020;
cplex_parameters["CPXERR_CANT_CLOSE_CHILD"] =      1021;
cplex_parameters["CPXERR_BAD_PROB_TYPE"] =         1022;
cplex_parameters["CPXERR_NOT_ONE_PROBLEM"] =       1023;
cplex_parameters["CPXERR_NOT_MILPCLASS"] =         1024;
cplex_parameters["CPXERR_STR_PARAM_TOO_LONG"] =    1026;
cplex_parameters["CPXERR_DECOMPRESSION"] =         1027;
cplex_parameters["CPXERR_BAD_PARAM_NAME"] =        1028;
cplex_parameters["CPXERR_NOT_MIQPCLASS"] =         1029;
cplex_parameters["CPXERR_NOT_FOR_QCP"] =           1031;

cplex_parameters["CPXERR_MSG_NO_CHANNEL"] =        1051;
cplex_parameters["CPXERR_MSG_NO_FILEPTR"] =        1052;
cplex_parameters["CPXERR_MSG_NO_FUNCTION"] =       1053;

cplex_parameters["CPXERR_PRESLV_INForUNBD"] =      1101;
cplex_parameters["CPXERR_PRESLV_NO_PROB"] =        1103;
cplex_parameters["CPXERR_PRESLV_ABORT"] =          1106;
cplex_parameters["CPXERR_PRESLV_BASIS_MEM"] =      1107;
cplex_parameters["CPXERR_PRESLV_COPYSOS"] =        1108;
cplex_parameters["CPXERR_PRESLV_COPYORDER"] =      1109;
cplex_parameters["CPXERR_PRESLV_SOLN_MIP"] =       1110;
cplex_parameters["CPXERR_PRESLV_SOLN_QP"] =        1111;
cplex_parameters["CPXERR_PRESLV_START_LP"] =       1112;
cplex_parameters["CPXERR_PRESLV_FAIL_BASIS"] =     1114;
cplex_parameters["CPXERR_PRESLV_NO_BASIS"] =       1115;
cplex_parameters["CPXERR_PRESLV_INF"] =            1117;
cplex_parameters["CPXERR_PRESLV_UNBD"] =           1118;
cplex_parameters["CPXERR_PRESLV_DUAL"] =           1119;
cplex_parameters["CPXERR_PRESLV_UNCRUSHFORM"] =    1120;
cplex_parameters["CPXERR_PRESLV_CRUSHFORM"] =      1121;
cplex_parameters["CPXERR_PRESLV_BAD_PARAM"] =      1122;
cplex_parameters["CPXERR_PRESLV_TIME_LIM"] =       1123;



/* Callable library miscellaneous routines */

cplex_parameters["CPXERR_INDEX_RANGE"] =           1200;
cplex_parameters["CPXERR_COL_INDEX_RANGE"] =       1201;
cplex_parameters["CPXERR_ROW_INDEX_RANGE"] =       1203;
cplex_parameters["CPXERR_INDEX_RANGE_LOW"] =       1205;
cplex_parameters["CPXERR_INDEX_RANGE_HIGH"] =      1206;
cplex_parameters["CPXERR_NEGATIVE_SURPLUS"] =      1207;
cplex_parameters["CPXERR_ARRAY_TOO_LONG"] =        1208;
cplex_parameters["CPXERR_NAME_CREATION"] =         1209;
cplex_parameters["CPXERR_NAME_NOT_FOUND"] =        1210;
cplex_parameters["CPXERR_NO_RHS_IN_OBJ"] =         1211;
cplex_parameters["CPXERR_BAD_SENSE"] =             1215;
cplex_parameters["CPXERR_NO_RNGVAL"] =             1216;
cplex_parameters["CPXERR_NO_SOLN"] =               1217;
cplex_parameters["CPXERR_NO_NAMES"] =              1219;
cplex_parameters["CPXERR_NOT_FIXED"] =             1221;
cplex_parameters["CPXERR_DUP_ENTRY"] =             1222;
cplex_parameters["CPXERR_NO_BARRIER_SOLN"] =       1223;
cplex_parameters["CPXERR_NULL_NAME"] =             1224;
cplex_parameters["CPXERR_NAN"] =                   1225;
cplex_parameters["CPXERR_ARRAY_NOT_ASCENDING"] =   1226;
cplex_parameters["CPXERR_COUNT_RANGE"] =           1227;
cplex_parameters["CPXERR_COUNT_OVERLAP"] =         1228;
cplex_parameters["CPXERR_BAD_LUB"] =               1229;
cplex_parameters["CPXERR_NODE_INDEX_RANGE"] =      1230;
cplex_parameters["CPXERR_ARC_INDEX_RANGE"] =       1231;
cplex_parameters["CPXERR_NO_DUAL_SOLN"] =          1232;
cplex_parameters["CPXERR_DBL_MAX"] =               1233;
cplex_parameters["CPXERR_THREAD_FAILED"] =         1234;


/* Simplex related */

cplex_parameters["CPXERR_INDEX_NOT_BASIC"] =       1251;
cplex_parameters["CPXERR_NEED_OPT_SOLN"] =         1252;
cplex_parameters["CPXERR_BAD_STATUS"] =            1253;
cplex_parameters["CPXERR_NOT_UNBOUNDED"] =         1254;
cplex_parameters["CPXERR_SBASE_INCOMPAT"] =        1255;
cplex_parameters["CPXERR_SINGULAR"] =              1256;
cplex_parameters["CPXERR_PRIIND"] =                1257;
cplex_parameters["CPXERR_NO_LU_FACTOR"] =          1258;
cplex_parameters["CPXERR_NO_SENSIT"] =             1260;
cplex_parameters["CPXERR_NO_BASIC_SOLN"] =         1261;
cplex_parameters["CPXERR_NO_BASIS"] =              1262;
cplex_parameters["CPXERR_ABORT_STRONGBRANCH"] =    1263;
cplex_parameters["CPXERR_NO_NORMS"] =              1264;
cplex_parameters["CPXERR_NOT_DUAL_UNBOUNDED"] =    1265;
cplex_parameters["CPXERR_TILIM_STRONGBRANCH"] =    1266;
cplex_parameters["CPXERR_BAD_PIVOT"] =             1267;
cplex_parameters["CPXERR_TILIM_CONDITION_NO"] =    1268;


/* Hybrid solvers */
cplex_parameters["CPXERR_BAD_METHOD"] =            1292;


/* For readers and writers */
cplex_parameters["CPXERR_NO_FILENAME"] =           1421;
cplex_parameters["CPXERR_FAIL_OPEN_WRITE"] =       1422;
cplex_parameters["CPXERR_FAIL_OPEN_READ"] =        1423;
cplex_parameters["CPXERR_BAD_FILETYPE"] =          1424;
cplex_parameters["CPXERR_XMLPARSE"] =              1425;


/* Common to LP, MPS, and related readers */

cplex_parameters["CPXERR_TOO_MANY_ROWS"] =         1431;
cplex_parameters["CPXERR_TOO_MANY_COLS"] =         1432;
cplex_parameters["CPXERR_TOO_MANY_COEFFS"] =       1433;
cplex_parameters["CPXERR_BAD_NUMBER"] =            1434;
cplex_parameters["CPXERR_BAD_EXPO_RANGE"] =        1435;
cplex_parameters["CPXERR_NO_OBJ_SENSE"] =          1436;
cplex_parameters["CPXERR_QCP_SENSE_FILE"] =        1437;
cplex_parameters["CPXERR_BAD_LAZY_UCUT"] =         1438;
cplex_parameters["CPXERR_BAD_INDCONSTR"] =         1439;

/* Common to MPS and related readers */

cplex_parameters["CPXERR_NO_NAME_SECTION"] =       1441;
cplex_parameters["CPXERR_BAD_SOS_TYPE"] =          1442;
cplex_parameters["CPXERR_COL_ROW_REPEATS"] =       1443;
cplex_parameters["CPXERR_RIM_ROW_REPEATS"] =       1444;
cplex_parameters["CPXERR_ROW_REPEATS"] =           1445;
cplex_parameters["CPXERR_COL_REPEATS"] =           1446;
cplex_parameters["CPXERR_RIM_REPEATS"] =           1447;
cplex_parameters["CPXERR_ROW_UNKNOWN"] =           1448;
cplex_parameters["CPXERR_COL_UNKNOWN"] =           1449;
cplex_parameters["CPXERR_NO_ROW_SENSE"] =          1453;
cplex_parameters["CPXERR_EXTRA_FX_BOUND"] =        1454;
cplex_parameters["CPXERR_EXTRA_FR_BOUND"] =        1455;
cplex_parameters["CPXERR_EXTRA_BV_BOUND"] =        1456;
cplex_parameters["CPXERR_BAD_BOUND_TYPE"] =        1457;
cplex_parameters["CPXERR_UP_BOUND_REPEATS"] =      1458;
cplex_parameters["CPXERR_LO_BOUND_REPEATS"] =      1459;
cplex_parameters["CPXERR_NO_BOUND_TYPE"] =         1460;
cplex_parameters["CPXERR_NO_QMATRIX_SECTION"] =    1461;
cplex_parameters["CPXERR_BAD_SECTION_ENDATA"] =    1462;
cplex_parameters["CPXERR_INT_TOO_BIG_INPUT"] =     1463;
cplex_parameters["CPXERR_NAME_TOO_LONG"] =         1464;
cplex_parameters["CPXERR_LINE_TOO_LONG"] =         1465;

/* Unique to MPS reader */

cplex_parameters["CPXERR_NO_ROWS_SECTION"] =       1471;
cplex_parameters["CPXERR_NO_COLUMNS_SECTION"] =    1472;
cplex_parameters["CPXERR_BAD_SECTION_BOUNDS"] =    1473;
cplex_parameters["CPXERR_RANGE_SECTION_ORDER"] =   1474;
cplex_parameters["CPXERR_BAD_SECTION_QMATRIX"] =   1475;
cplex_parameters["CPXERR_NO_OBJECTIVE"] =          1476;
cplex_parameters["CPXERR_ROW_REPEAT_PRINT"] =      1477;
cplex_parameters["CPXERR_COL_REPEAT_PRINT"] =      1478;
cplex_parameters["CPXERR_RIMNZ_REPEATS"] =         1479;
cplex_parameters["CPXERR_EXTRA_INTORG"] =          1480;
cplex_parameters["CPXERR_EXTRA_INTEND"] =          1481;
cplex_parameters["CPXERR_EXTRA_SOSORG"] =          1482;
cplex_parameters["CPXERR_EXTRA_SOSEND"] =          1483;
cplex_parameters["CPXERR_TOO_MANY_RIMS"] =         1484;
cplex_parameters["CPXERR_TOO_MANY_RIMNZ"] =        1485;
cplex_parameters["CPXERR_NO_ROW_NAME"] =           1486;
cplex_parameters["CPXERR_BAD_OBJ_SENSE"] =         1487;

/* BAS files */

cplex_parameters["CPXERR_BAS_FILE_SHORT"] =        1550;
cplex_parameters["CPXERR_BAD_INDICATOR"] =         1551;
cplex_parameters["CPXERR_NO_ENDATA"] =             1552;
cplex_parameters["CPXERR_FILE_ENTRIES"] =          1553;
cplex_parameters["CPXERR_SBASE_ILLEGAL"] =         1554;
cplex_parameters["CPXERR_BAS_FILE_SIZE"] =         1555;
cplex_parameters["CPXERR_NO_VECTOR_SOLN"] =        1556;

/* SAV files */

cplex_parameters["CPXERR_NOT_SAV_FILE"] =          1560;
cplex_parameters["CPXERR_SAV_FILE_DATA"] =         1561;
cplex_parameters["CPXERR_SAV_FILE_WRITE"] =        1562;
cplex_parameters["CPXERR_FILE_FORMAT"] =           1563;

/* LP reader errors */

cplex_parameters["CPXERR_ADJ_SIGNS"] =             1602;
cplex_parameters["CPXERR_RHS_IN_OBJ"] =            1603;
cplex_parameters["CPXERR_ADJ_SIGN_SENSE"] =        1604;
cplex_parameters["CPXERR_QUAD_IN_ROW"] =           1605;
cplex_parameters["CPXERR_ADJ_SIGN_QUAD"] =         1606;
cplex_parameters["CPXERR_NO_OPERATOR"] =           1607;
cplex_parameters["CPXERR_NO_OP_OR_SENSE"] =        1608;
cplex_parameters["CPXERR_NO_ID_FIRST"] =           1609;
cplex_parameters["CPXERR_NO_RHS_COEFF"] =          1610;
cplex_parameters["CPXERR_NO_NUMBER_FIRST"] =       1611;
cplex_parameters["CPXERR_NO_QUAD_EXP"] =           1612;
cplex_parameters["CPXERR_QUAD_EXP_NOT_2"] =        1613;
cplex_parameters["CPXERR_NO_QP_OPERATOR"] =        1614;
cplex_parameters["CPXERR_NO_NUMBER"] =             1615;
cplex_parameters["CPXERR_NO_ID"] =                 1616;
cplex_parameters["CPXERR_BAD_ID"] =                1617;
cplex_parameters["CPXERR_BAD_EXPONENT"] =          1618;
cplex_parameters["CPXERR_Q_DIVISOR"] =             1619;
cplex_parameters["CPXERR_NO_BOUND_SENSE"] =        1621;
cplex_parameters["CPXERR_BAD_BOUND_SENSE"] =       1622;
cplex_parameters["CPXERR_NO_NUMBER_BOUND"] =       1623;
cplex_parameters["CPXERR_NO_SOS_SEPARATOR"] =      1627;
cplex_parameters["CPXERR_INVALID_NUMBER"] =        1650;

/* Parameter file read errors */

cplex_parameters["CPXERR_PRM_DATA"] =              1660;
cplex_parameters["CPXERR_PRM_HEADER"] =            1661;


/* Conflict error codes */

cplex_parameters["CPXERR_NO_CONFLICT"] =          1719;

cplex_parameters["CPXERR_CONFLICT_UNSTABLE"] =    1720;



/* Temporary file errors */

cplex_parameters["CPXERR_WORK_FILE_OPEN"] =          1801;
cplex_parameters["CPXERR_WORK_FILE_READ"] =          1802;
cplex_parameters["CPXERR_WORK_FILE_WRITE"] =         1803;
cplex_parameters["CPXERR_IN_INFOCALLBACK"] =         1804;
cplex_parameters["CPXERR_MIPSEARCH_WITH_CALLBACKS"] = 1805;
cplex_parameters["CPXERR_LP_NOT_IN_ENVIRONMENT"] =1806;

cplex_parameters["CPXERR_PARAM_INCOMPATIBLE"] =        1807;

/* Private errors are 2000-2999 (cpxpriv.h)
   MIP errors are 3000-3999 (mipdefs.h)
   Barrier errors are 4000-4999 (bardefs.h)
   QP errors are 5000-5999 (qpdefs.h) */

/* Licensor errors */

cplex_parameters["CPXERR_LICENSE_MIN"] =          32000;
cplex_parameters["CPXERR_ILOG_LICENSE"] =         32201;

cplex_parameters["CPXERR_NO_MIP_LIC"] =           32301;
cplex_parameters["CPXERR_NO_BARRIER_LIC"] =       32302;
cplex_parameters["CPXERR_NO_MIQP_LIC"] =          32305;
cplex_parameters["CPXERR_BADLDWID"] =             32018;

cplex_parameters["CPXERR_BADPRODUCT"] =           32023;

cplex_parameters["CPXERR_ALGNOTLICENSED"] =       32024;
cplex_parameters["CPXERR_LICENSE_MAX"] =          32999;

/* Generic constants */

cplex_parameters["CPX_ON"] =                          1;
cplex_parameters["CPX_OFF"] =                         0;
cplex_parameters["CPX_MAX"] =                        -1;
cplex_parameters["CPX_MIN"] =                         1;

/* Pricing options */

cplex_parameters["CPX_PPRIIND_PARTIAL"] =            -1;
cplex_parameters["CPX_PPRIIND_AUTO"] =                0;
cplex_parameters["CPX_PPRIIND_DEVEX"] =               1;
cplex_parameters["CPX_PPRIIND_STEEP"] =               2;
cplex_parameters["CPX_PPRIIND_STEEPQSTART"] =         3;
cplex_parameters["CPX_PPRIIND_FULL"] =                4;

cplex_parameters["CPX_DPRIIND_AUTO"] =                0;
cplex_parameters["CPX_DPRIIND_FULL"] =                1;
cplex_parameters["CPX_DPRIIND_STEEP"] =               2;
cplex_parameters["CPX_DPRIIND_FULLSTEEP"] =           3;
cplex_parameters["CPX_DPRIIND_STEEPQSTART"] =         4;
cplex_parameters["CPX_DPRIIND_DEVEX"] =               5;

/* PARALLELMODE values  */

cplex_parameters["CPX_PARALLEL_DETERMINISTIC"] =     1;
cplex_parameters["CPX_PARALLEL_AUTO"] =              0;
cplex_parameters["CPX_PARALLEL_OPPORTUNISTIC"] =    -1;

/* Values for CPX_PARAM_WRITELEVEL */

cplex_parameters["CPX_WRITELEVEL_AUTO"] =                0;
cplex_parameters["CPX_WRITELEVEL_ALLVARS"] =             1;
cplex_parameters["CPX_WRITELEVEL_DISCRETEVARS"] =        2;
cplex_parameters["CPX_WRITELEVEL_NONZEROVARS"] =         3;
cplex_parameters["CPX_WRITELEVEL_NONZERODISCRETEVARS"] = 4 ;

/* LP/QP solution algorithms, used as possible values for
   CPX_PARAM_LPMETHOD/CPX_PARAM_QPMETHOD/CPX_PARAM_BARCROSSALG/
   CPXgetmethod/... */

cplex_parameters["CPX_ALG_NONE"] =                    -1;
cplex_parameters["CPX_ALG_AUTOMATIC"] =                0;
cplex_parameters["CPX_ALG_PRIMAL"] =                   1;
cplex_parameters["CPX_ALG_DUAL"] =                     2;
cplex_parameters["CPX_ALG_NET"] =                      3;
cplex_parameters["CPX_ALG_BARRIER"] =                  4;
cplex_parameters["CPX_ALG_SIFTING"] =                  5;
cplex_parameters["CPX_ALG_CONCURRENT"] =               6;
cplex_parameters["CPX_ALG_BAROPT"] =                   7;
cplex_parameters["CPX_ALG_PIVOTIN"] =                  8;
cplex_parameters["CPX_ALG_PIVOTOUT"] =                 9;
cplex_parameters["CPX_ALG_PIVOT"] =                   10;
cplex_parameters["CPX_ALG_FEASOPT"] =                 11;
cplex_parameters["CPX_ALG_MIP"] =                     12;
cplex_parameters["CPX_ALG_ROBUST"] =                  13;

/* Basis status values */

cplex_parameters["CPX_AT_LOWER"] =                    0;
cplex_parameters["CPX_BASIC"] =                       1;
cplex_parameters["CPX_AT_UPPER"] =                    2;
cplex_parameters["CPX_FREE_SUPER"] =                  3;

/* Used in pivoting interface */

cplex_parameters["CPX_NO_VARIABLE"] =        2100000000;


/* Variable types for ctype array */

cplex_parameters["CPX_CONTINUOUS"] =                  'C';
cplex_parameters["CPX_BINARY"] =                      'B';
cplex_parameters["CPX_INTEGER"] =                     'I';
cplex_parameters["CPX_SEMICONT"] =                    'S';
cplex_parameters["CPX_SEMIINT"] =                     'N';

/* PREREDUCE settings */

cplex_parameters["CPX_PREREDUCE_PRIMALANDDUAL"] =      3;
cplex_parameters["CPX_PREREDUCE_DUALONLY"] =           2;
cplex_parameters["CPX_PREREDUCE_PRIMALONLY"] =         1;
cplex_parameters["CPX_PREREDUCE_NOPRIMALORDUAL"] =     0;

/* Conflict statuses */

cplex_parameters["CPX_STAT_CONFLICT_FEASIBLE"] =                  30;

cplex_parameters["CPX_STAT_CONFLICT_MINIMAL"] =                   31;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_CONTRADICTION"] =       32;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_TIME_LIM"] =            33;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_IT_LIM"] =              34;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_NODE_LIM"] =            35;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_OBJ_LIM"] =             36;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_MEM_LIM"] =             37;

cplex_parameters["CPX_STAT_CONFLICT_ABORT_USER"] =                38;


/* Conflict status values */

cplex_parameters["CPX_CONFLICT_EXCLUDED"] =         -1;
cplex_parameters["CPX_CONFLICT_POSSIBLE_MEMBER"] =   0;
cplex_parameters["CPX_CONFLICT_POSSIBLE_LB"] =       1;
cplex_parameters["CPX_CONFLICT_POSSIBLE_UB"] =       2;
cplex_parameters["CPX_CONFLICT_MEMBER"] =            3;
cplex_parameters["CPX_CONFLICT_LB"] =                4;
cplex_parameters["CPX_CONFLICT_UB"] =                5;

/* Problem Types
   Types 4, 9, and 12 are internal, the others are for users */

cplex_parameters["CPXPROB_LP"] =                      0;
cplex_parameters["CPXPROB_MILP"] =                    1;
cplex_parameters["CPXPROB_FIXEDMILP"] =               3;
cplex_parameters["CPXPROB_NODELP"] =                  4;
cplex_parameters["CPXPROB_QP"] =                      5;
cplex_parameters["CPXPROB_MIQP"] =                    7;
cplex_parameters["CPXPROB_FIXEDMIQP"] =               8;
cplex_parameters["CPXPROB_NODEQP"] =                  9;
cplex_parameters["CPXPROB_QCP"] =                    10;
cplex_parameters["CPXPROB_MIQCP"] =                  11;
cplex_parameters["CPXPROB_NODEQCP"] =                12;


/* CPLEX Parameter numbers */

cplex_parameters["CPX_PARAM_ADVIND"] =              1001;
cplex_parameters["CPX_PARAM_AGGFILL"] =             1002;
cplex_parameters["CPX_PARAM_AGGIND"] =              1003;
cplex_parameters["CPX_PARAM_BASINTERVAL"] =         1004;
cplex_parameters["CPX_PARAM_CFILEMUL"] =            1005;
cplex_parameters["CPX_PARAM_CLOCKTYPE"] =           1006;
cplex_parameters["CPX_PARAM_CRAIND"] =              1007;
cplex_parameters["CPX_PARAM_DEPIND"] =              1008;
cplex_parameters["CPX_PARAM_DPRIIND"] =             1009;
cplex_parameters["CPX_PARAM_PRICELIM"] =            1010;
cplex_parameters["CPX_PARAM_EPMRK"] =               1013;
cplex_parameters["CPX_PARAM_EPOPT"] =               1014;
cplex_parameters["CPX_PARAM_EPPER"] =               1015;
cplex_parameters["CPX_PARAM_EPRHS"] =               1016;
cplex_parameters["CPX_PARAM_FASTMIP"] =             1017;
cplex_parameters["CPX_PARAM_SIMDISPLAY"] =          1019;
cplex_parameters["CPX_PARAM_ITLIM"] =               1020;
cplex_parameters["CPX_PARAM_ROWREADLIM"] =          1021;
cplex_parameters["CPX_PARAM_NETFIND"] =             1022;
cplex_parameters["CPX_PARAM_COLREADLIM"] =          1023;
cplex_parameters["CPX_PARAM_NZREADLIM"] =           1024;
cplex_parameters["CPX_PARAM_OBJLLIM"] =             1025;
cplex_parameters["CPX_PARAM_OBJULIM"] =             1026;
cplex_parameters["CPX_PARAM_PERIND"] =              1027;
cplex_parameters["CPX_PARAM_PERLIM"] =              1028;
cplex_parameters["CPX_PARAM_PPRIIND"] =             1029;
cplex_parameters["CPX_PARAM_PREIND"] =              1030;
cplex_parameters["CPX_PARAM_REINV"] =               1031;
cplex_parameters["CPX_PARAM_REVERSEIND"] =          1032;
cplex_parameters["CPX_PARAM_RFILEMUL"] =            1033;
cplex_parameters["CPX_PARAM_SCAIND"] =              1034;
cplex_parameters["CPX_PARAM_SCRIND"] =              1035;
cplex_parameters["CPX_PARAM_SINGLIM"] =             1037;
cplex_parameters["CPX_PARAM_SINGTOL"] =             1038;
cplex_parameters["CPX_PARAM_TILIM"] =               1039;
cplex_parameters["CPX_PARAM_XXXIND"] =              1041;
cplex_parameters["CPX_PARAM_PREDUAL"] =             1044;
cplex_parameters["CPX_PARAM_EPOPT_H"] =             1049;
cplex_parameters["CPX_PARAM_EPRHS_H"] =             1050;
cplex_parameters["CPX_PARAM_PREPASS"] =             1052;
cplex_parameters["CPX_PARAM_DATACHECK"] =           1056;
cplex_parameters["CPX_PARAM_REDUCE"] =              1057;
cplex_parameters["CPX_PARAM_PRELINEAR"] =           1058;
cplex_parameters["CPX_PARAM_LPMETHOD"] =            1062;
cplex_parameters["CPX_PARAM_QPMETHOD"] =            1063;
cplex_parameters["CPX_PARAM_WORKDIR"] =             1064;
cplex_parameters["CPX_PARAM_WORKMEM"] =             1065;
cplex_parameters["CPX_PARAM_THREADS"] =             1067;
cplex_parameters["CPX_PARAM_CONFLICTDISPLAY"] =     1074;
cplex_parameters["CPX_PARAM_SIFTDISPLAY"] =         1076;
cplex_parameters["CPX_PARAM_SIFTALG"] =             1077;
cplex_parameters["CPX_PARAM_SIFTITLIM"] =           1078;
cplex_parameters["CPX_PARAM_MPSLONGNUM"] =          1081;
cplex_parameters["CPX_PARAM_MEMORYEMPHASIS"] =      1082;
cplex_parameters["CPX_PARAM_NUMERICALEMPHASIS"] =   1083;
cplex_parameters["CPX_PARAM_FEASOPTMODE"] =         1084;
cplex_parameters["CPX_PARAM_PARALLELMODE"] =        1109;
cplex_parameters["CPX_PARAM_TUNINGMEASURE"] =       1110;
cplex_parameters["CPX_PARAM_TUNINGREPEAT"] =        1111;
cplex_parameters["CPX_PARAM_TUNINGTILIM"] =         1112;
cplex_parameters["CPX_PARAM_TUNINGDISPLAY"] =       1113;
cplex_parameters["CPX_PARAM_WRITELEVEL"] =          1114;

/* Barrier is in bardefs.h, MIP is in mipdefs.h, QP is in qpdefs.h */

cplex_parameters["CPX_PARAM_ALL_MIN"] =             1000;
cplex_parameters["CPX_PARAM_ALL_MAX"] =             6000;

/* Callback values for wherefrom   */

cplex_parameters["CPX_CALLBACK_PRIMAL"] =              1;
cplex_parameters["CPX_CALLBACK_DUAL"] =                2;
cplex_parameters["CPX_CALLBACK_NETWORK"] =             3;
cplex_parameters["CPX_CALLBACK_PRIMAL_CROSSOVER"] =    4;
cplex_parameters["CPX_CALLBACK_DUAL_CROSSOVER"] =      5;
cplex_parameters["CPX_CALLBACK_BARRIER"] =             6;
cplex_parameters["CPX_CALLBACK_PRESOLVE"] =            7;
cplex_parameters["CPX_CALLBACK_QPBARRIER"] =           8;
cplex_parameters["CPX_CALLBACK_QPSIMPLEX"] =           9;
cplex_parameters["CPX_CALLBACK_TUNING"] =             10;
/* Be sure to check the MIP values */

/* Values for getcallbackinfo function */

cplex_parameters["CPX_CALLBACK_INFO_PRIMAL_OBJ"] =           1;
cplex_parameters["CPX_CALLBACK_INFO_DUAL_OBJ"] =             2;
cplex_parameters["CPX_CALLBACK_INFO_PRIMAL_INFMEAS"] =       3;
cplex_parameters["CPX_CALLBACK_INFO_DUAL_INFMEAS"] =         4;
cplex_parameters["CPX_CALLBACK_INFO_PRIMAL_FEAS"] =          5;
cplex_parameters["CPX_CALLBACK_INFO_DUAL_FEAS"] =            6;
cplex_parameters["CPX_CALLBACK_INFO_ITCOUNT"] =              7;
cplex_parameters["CPX_CALLBACK_INFO_CROSSOVER_PPUSH"] =      8;
cplex_parameters["CPX_CALLBACK_INFO_CROSSOVER_PEXCH"] =      9;
cplex_parameters["CPX_CALLBACK_INFO_CROSSOVER_DPUSH"] =     10;
cplex_parameters["CPX_CALLBACK_INFO_CROSSOVER_DEXCH"] =     11;
cplex_parameters["CPX_CALLBACK_INFO_CROSSOVER_SBCNT"] =     12;
cplex_parameters["CPX_CALLBACK_INFO_PRESOLVE_ROWSGONE"] =   13;
cplex_parameters["CPX_CALLBACK_INFO_PRESOLVE_COLSGONE"] =   14;
cplex_parameters["CPX_CALLBACK_INFO_PRESOLVE_AGGSUBST"] =   15;
cplex_parameters["CPX_CALLBACK_INFO_PRESOLVE_COEFFS"] =     16;
cplex_parameters["CPX_CALLBACK_INFO_USER_PROBLEM"] =        17;
cplex_parameters["CPX_CALLBACK_INFO_TUNING_PROGRESS"] =     18;
cplex_parameters["CPX_CALLBACK_INFO_ENDTIME"] =             19;
/* Be sure to check the MIP values */


/* Values for CPX_PARAM_TUNINGMEASURE */
cplex_parameters["CPX_TUNE_AVERAGE"] = 1;
cplex_parameters["CPX_TUNE_MINMAX"] =  2;


/* Values for incomplete tuning */
cplex_parameters["CPX_TUNE_ABORT"] =    1;
cplex_parameters["CPX_TUNE_TILIM"] =    2;


/* Quality query identifiers */

cplex_parameters["CPX_MAX_PRIMAL_INFEAS"] =         1;
cplex_parameters["CPX_MAX_SCALED_PRIMAL_INFEAS"] =  2;
cplex_parameters["CPX_SUM_PRIMAL_INFEAS"] =         3;
cplex_parameters["CPX_SUM_SCALED_PRIMAL_INFEAS"] =  4;
cplex_parameters["CPX_MAX_DUAL_INFEAS"] =           5;
cplex_parameters["CPX_MAX_SCALED_DUAL_INFEAS"] =    6;
cplex_parameters["CPX_SUM_DUAL_INFEAS"] =           7;
cplex_parameters["CPX_SUM_SCALED_DUAL_INFEAS"] =    8;
cplex_parameters["CPX_MAX_INT_INFEAS"] =            9;
cplex_parameters["CPX_SUM_INT_INFEAS"] =            10;
cplex_parameters["CPX_MAX_PRIMAL_RESIDUAL"] =       11;
cplex_parameters["CPX_MAX_SCALED_PRIMAL_RESIDUAL"] =12;
cplex_parameters["CPX_SUM_PRIMAL_RESIDUAL"] =       13;
cplex_parameters["CPX_SUM_SCALED_PRIMAL_RESIDUAL"] =14;
cplex_parameters["CPX_MAX_DUAL_RESIDUAL"] =         15;
cplex_parameters["CPX_MAX_SCALED_DUAL_RESIDUAL"] =  16;
cplex_parameters["CPX_SUM_DUAL_RESIDUAL"] =         17;
cplex_parameters["CPX_SUM_SCALED_DUAL_RESIDUAL"] =  18;
cplex_parameters["CPX_MAX_COMP_SLACK"] =            19;
cplex_parameters["CPX_SUM_COMP_SLACK"] =            21;
cplex_parameters["CPX_MAX_X"] =                     23;
cplex_parameters["CPX_MAX_SCALED_X"] =              24;
cplex_parameters["CPX_MAX_PI"] =                    25;
cplex_parameters["CPX_MAX_SCALED_PI"] =             26;
cplex_parameters["CPX_MAX_SLACK"] =                 27;
cplex_parameters["CPX_MAX_SCALED_SLACK"] =          28;
cplex_parameters["CPX_MAX_RED_COST"] =              29;
cplex_parameters["CPX_MAX_SCALED_RED_COST"] =       30;
cplex_parameters["CPX_SUM_X"] =                     31;
cplex_parameters["CPX_SUM_SCALED_X"] =              32;
cplex_parameters["CPX_SUM_PI"] =                    33;
cplex_parameters["CPX_SUM_SCALED_PI"] =             34;
cplex_parameters["CPX_SUM_SLACK"] =                 35;
cplex_parameters["CPX_SUM_SCALED_SLACK"] =          36;
cplex_parameters["CPX_SUM_RED_COST"] =              37;
cplex_parameters["CPX_SUM_SCALED_RED_COST"] =       38;
cplex_parameters["CPX_KAPPA"] =                     39;
cplex_parameters["CPX_OBJ_GAP"] =                   40;
cplex_parameters["CPX_DUAL_OBJ"] =                  41;
cplex_parameters["CPX_PRIMAL_OBJ"] =                42;
cplex_parameters["CPX_MAX_QCPRIMAL_RESIDUAL"] =     43;
cplex_parameters["CPX_SUM_QCPRIMAL_RESIDUAL"] =     44;
cplex_parameters["CPX_MAX_QCSLACK_INFEAS"] =        45;
cplex_parameters["CPX_SUM_QCSLACK_INFEAS"] =        46;
cplex_parameters["CPX_MAX_QCSLACK"] =               47;
cplex_parameters["CPX_SUM_QCSLACK"] =               48;
cplex_parameters["CPX_MAX_INDSLACK_INFEAS"] =       49;
cplex_parameters["CPX_SUM_INDSLACK_INFEAS"] =       50;
cplex_parameters["CPX_EXACT_KAPPA"] =               51;


/* IIS errors */

cplex_parameters["CPXERR_IIS_NO_INFO"] =           1701;
cplex_parameters["CPXERR_IIS_NO_SOLN"] =           1702;
cplex_parameters["CPXERR_IIS_FEAS"] =              1703;
cplex_parameters["CPXERR_IIS_NOT_INFEAS"] =        1704;
cplex_parameters["CPXERR_IIS_OPT_INFEAS"] =        1705;
cplex_parameters["CPXERR_IIS_DEFAULT"] =           1706;
cplex_parameters["CPXERR_IIS_NO_BASIC"] =          1707;
cplex_parameters["CPXERR_IIS_NO_LOAD"] =           1709;
cplex_parameters["CPXERR_IIS_SUB_OBJ_LIM"] =       1710;
cplex_parameters["CPXERR_IIS_SUB_IT_LIM"] =        1711;
cplex_parameters["CPXERR_IIS_SUB_TIME_LIM"] =      1712;
cplex_parameters["CPXERR_IIS_NUM_BEST"] =          1713;
cplex_parameters["CPXERR_IIS_SUB_ABORT"] =         1714;

/* Infeasibility Finder return values */

cplex_parameters["CPXIIS_COMPLETE"] =                 1;
cplex_parameters["CPXIIS_PARTIAL"] =                  2;

/* Infeasibility Finder display values */

cplex_parameters["CPXIIS_TERSE"] =                     1;
cplex_parameters["CPXIIS_VERBOSE"] =                   2;

/* Infeasibility Finder row and column statuses */

cplex_parameters["CPXIIS_AT_LOWER"] =                  0;
cplex_parameters["CPXIIS_FIXED"] =                     1;
cplex_parameters["CPXIIS_AT_UPPER"] =                  2;

/* --------------------------------------------------------------------------
 * File: bardefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 1994, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */


cplex_parameters["CPX_STAT_OPTIMAL_FACE_UNBOUNDED"] =  20;
cplex_parameters["CPX_STAT_ABORT_PRIM_OBJ_LIM"] =      21;
cplex_parameters["CPX_STAT_ABORT_DUAL_OBJ_LIM"] =      22;


/* Barrier parameters */

cplex_parameters["CPX_PARAM_BARDSTART"] =           3001;
cplex_parameters["CPX_PARAM_BAREPCOMP"] =           3002;
cplex_parameters["CPX_PARAM_BARGROWTH"] =           3003;
cplex_parameters["CPX_PARAM_BAROBJRNG"] =           3004;
cplex_parameters["CPX_PARAM_BARPSTART"] =           3005;
cplex_parameters["CPX_PARAM_BARALG"] =              3007;
cplex_parameters["CPX_PARAM_BARCOLNZ"] =            3009;
cplex_parameters["CPX_PARAM_BARDISPLAY"] =          3010;
cplex_parameters["CPX_PARAM_BARITLIM"] =            3012;
cplex_parameters["CPX_PARAM_BARMAXCOR"] =           3013;
cplex_parameters["CPX_PARAM_BARORDER"] =            3014;
cplex_parameters["CPX_PARAM_BARSTARTALG"] =         3017;
cplex_parameters["CPX_PARAM_BARCROSSALG"] =         3018;
cplex_parameters["CPX_PARAM_BARQCPEPCOMP"] =        3020;

/* Optimizing Problems */

cplex_parameters["CPX_BARORDER_AUTO"] =0;
cplex_parameters["CPX_BARORDER_AMD"] = 1;
cplex_parameters["CPX_BARORDER_AMF"] = 2;
cplex_parameters["CPX_BARORDER_ND"] =  3;



/* --------------------------------------------------------------------------
 * File: mipdefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 1991, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

/* MIP error codes */

cplex_parameters["CPXERR_NOT_MIP"] =               3003;
cplex_parameters["CPXERR_BAD_PRIORITY"] =          3006;
cplex_parameters["CPXERR_ORDER_BAD_DIRECTION"] =   3007;
cplex_parameters["CPXERR_ARRAY_BAD_SOS_TYPE"] =    3009;
cplex_parameters["CPXERR_UNIQUE_WEIGHTS"] =        3010;
cplex_parameters["CPXERR_BAD_DIRECTION"] =         3012;
cplex_parameters["CPXERR_NO_SOS"] =                3015;
cplex_parameters["CPXERR_NO_ORDER"] =              3016;
cplex_parameters["CPXERR_INT_TOO_BIG"] =           3018;
cplex_parameters["CPXERR_SUBPROB_SOLVE"] =         3019;
cplex_parameters["CPXERR_NO_MIPSTART"] =           3020;
cplex_parameters["CPXERR_BAD_CTYPE"] =             3021;
cplex_parameters["CPXERR_NO_INT_X"] =              3023;

cplex_parameters["CPXERR_NO_SOLNPOOL"] =           3024;

cplex_parameters["CPXERR_MISS_SOS_TYPE"] =         3301;
       
cplex_parameters["CPXERR_NO_TREE"] =               3412;
cplex_parameters["CPXERR_TREE_MEMORY_LIMIT"] =     3413;
cplex_parameters["CPXERR_FILTER_VARIABLE_TYPE"] =  3414;
       
cplex_parameters["CPXERR_NODE_ON_DISK"] =          3504;
       
cplex_parameters["CPXERR_PTHREAD_MUTEX_INIT"] =    3601;
cplex_parameters["CPXERR_PTHREAD_CREATE"] =        3603;

/* MIP emphasis settings */

cplex_parameters["CPX_MIPEMPHASIS_BALANCED"] =    0;
cplex_parameters["CPX_MIPEMPHASIS_FEASIBILITY"] = 1;
cplex_parameters["CPX_MIPEMPHASIS_OPTIMALITY"] =  2;
cplex_parameters["CPX_MIPEMPHASIS_BESTBOUND"] =   3;
cplex_parameters["CPX_MIPEMPHASIS_HIDDENFEAS"] =  4;

/* Values for sostype and branch type */

cplex_parameters["CPX_TYPE_VAR"] =                   '0';
cplex_parameters["CPX_TYPE_SOS1"] =                  '1';
cplex_parameters["CPX_TYPE_SOS2"] =                  '2';
cplex_parameters["CPX_TYPE_USER"] =                  'X';
cplex_parameters["CPX_TYPE_ANY"] =                   'A';

/* Variable selection values */

cplex_parameters["CPX_VARSEL_MININFEAS"] =           -1;
cplex_parameters["CPX_VARSEL_DEFAULT"] =              0;
cplex_parameters["CPX_VARSEL_MAXINFEAS"] =            1;
cplex_parameters["CPX_VARSEL_PSEUDO"] =               2;
cplex_parameters["CPX_VARSEL_STRONG"] =               3;
cplex_parameters["CPX_VARSEL_PSEUDOREDUCED"] =        4;

/* Node selection values */

cplex_parameters["CPX_NODESEL_DFS"] =                 0;
cplex_parameters["CPX_NODESEL_BESTBOUND"] =           1;
cplex_parameters["CPX_NODESEL_BESTEST"] =             2;
cplex_parameters["CPX_NODESEL_BESTEST_ALT"] =         3;

/* Values for generated priority order */

cplex_parameters["CPX_MIPORDER_COST"] =               1 ;
cplex_parameters["CPX_MIPORDER_BOUNDS"] =             2 ;
cplex_parameters["CPX_MIPORDER_SCALEDCOST"] =         3;

/* Values for direction array */

cplex_parameters["CPX_BRANCH_GLOBAL"] =               0;
cplex_parameters["CPX_BRANCH_DOWN"] =                -1;
cplex_parameters["CPX_BRANCH_UP"] =                   1;

/* Values for CPX_PARAM_BRDIR */

cplex_parameters["CPX_BRDIR_DOWN"] =                 -1;
cplex_parameters["CPX_BRDIR_AUTO"] =                  0;
cplex_parameters["CPX_BRDIR_UP"] =                    1;
 
/* Values for cuttype in CPXgetnumcuts */
/* !!! Update mipnodelog.c::cutname[] if you change any of these !!! */

cplex_parameters["CPX_CUT_COVER"] =     0;
cplex_parameters["CPX_CUT_GUBCOVER"] =  1;
cplex_parameters["CPX_CUT_FLOWCOVER"] = 2;
cplex_parameters["CPX_CUT_CLIQUE"] =    3;
cplex_parameters["CPX_CUT_FRAC"] =      4;
cplex_parameters["CPX_CUT_MIR"] =       5;
cplex_parameters["CPX_CUT_FLOWPATH"] =  6;
cplex_parameters["CPX_CUT_DISJ"] =      7;
cplex_parameters["CPX_CUT_IMPLBD"] =    8;
cplex_parameters["CPX_CUT_ZEROHALF"] =  9;
cplex_parameters["CPX_CUT_MCF"] =       10;
cplex_parameters["CPX_CUT_LOCALCOVER"] =11;
cplex_parameters["CPX_CUT_TIGHTEN"] =   12;
cplex_parameters["CPX_CUT_OBJDISJ"] =   13;
cplex_parameters["CPX_CUT_USER"] =      14;
cplex_parameters["CPX_CUT_TABLE"] =     15;
cplex_parameters["CPX_CUT_SOLNPOOL"] =  16;
cplex_parameters["CPX_CUT_NUM_TYPES"] = 17;


/* Values for CPX_PARAM_MIPSEARCH */

cplex_parameters["CPX_MIPSEARCH_AUTO"] =              0;
cplex_parameters["CPX_MIPSEARCH_TRADITIONAL"] =       1;
cplex_parameters["CPX_MIPSEARCH_DYNAMIC"] =           2;
 
/* Effort levels for MIP starts */

cplex_parameters["CPX_MIPSTART_AUTO"] =        0;
cplex_parameters["CPX_MIPSTART_CHECKFEAS"] =   1;
cplex_parameters["CPX_MIPSTART_SOLVEFIXED"] =  2;
cplex_parameters["CPX_MIPSTART_SOLVEMIP"] =    3;
cplex_parameters["CPX_MIPSTART_REPAIR"] =      4;


/* MIP Problem status codes */

cplex_parameters["CPXMIP_OPTIMAL"] =                101;
cplex_parameters["CPXMIP_OPTIMAL_TOL"] =            102;
cplex_parameters["CPXMIP_INFEASIBLE"] =             103;
cplex_parameters["CPXMIP_SOL_LIM"] =                104;
cplex_parameters["CPXMIP_NODE_LIM_FEAS"] =          105;
cplex_parameters["CPXMIP_NODE_LIM_INFEAS"] =        106;
cplex_parameters["CPXMIP_TIME_LIM_FEAS"] =          107;
cplex_parameters["CPXMIP_TIME_LIM_INFEAS"] =        108;
cplex_parameters["CPXMIP_FAIL_FEAS"] =              109;
cplex_parameters["CPXMIP_FAIL_INFEAS"] =            110;
cplex_parameters["CPXMIP_MEM_LIM_FEAS"] =           111;
cplex_parameters["CPXMIP_MEM_LIM_INFEAS"] =         112;
cplex_parameters["CPXMIP_ABORT_FEAS"] =             113;
cplex_parameters["CPXMIP_ABORT_INFEAS"] =           114;
cplex_parameters["CPXMIP_OPTIMAL_INFEAS"] =         115;
cplex_parameters["CPXMIP_FAIL_FEAS_NO_TREE"] =      116;
cplex_parameters["CPXMIP_FAIL_INFEAS_NO_TREE"] =    117;
cplex_parameters["CPXMIP_UNBOUNDED"] =              118;
cplex_parameters["CPXMIP_INForUNBD"] =              119;

cplex_parameters["CPXMIP_FEASIBLE_RELAXED_SUM"] =   120;

cplex_parameters["CPXMIP_OPTIMAL_RELAXED_SUM"] =    121;

cplex_parameters["CPXMIP_FEASIBLE_RELAXED_INF"] =   122;

cplex_parameters["CPXMIP_OPTIMAL_RELAXED_INF"] =    123;

cplex_parameters["CPXMIP_FEASIBLE_RELAXED_QUAD"] =  124;

cplex_parameters["CPXMIP_OPTIMAL_RELAXED_QUAD"] =   125;

cplex_parameters["CPXMIP_ABORT_RELAXED"] =   126;


cplex_parameters["CPXMIP_FEASIBLE"] =127;

cplex_parameters["CPXMIP_POPULATESOL_LIM"] =128;

cplex_parameters["CPXMIP_OPTIMAL_POPULATED"] =129;

cplex_parameters["CPXMIP_OPTIMAL_POPULATED_TOL"] =130;


/* Callback values for wherefrom */

cplex_parameters["CPX_CALLBACK_MIP"] =              101;
cplex_parameters["CPX_CALLBACK_MIP_BRANCH"] =       102;
cplex_parameters["CPX_CALLBACK_MIP_NODE"] =         103;
cplex_parameters["CPX_CALLBACK_MIP_HEURISTIC"] =    104;
cplex_parameters["CPX_CALLBACK_MIP_SOLVE"] =        105;
cplex_parameters["CPX_CALLBACK_MIP_CUT"] =          106;
cplex_parameters["CPX_CALLBACK_MIP_PROBE"] =        107;
cplex_parameters["CPX_CALLBACK_MIP_FRACCUT"] =      108;
cplex_parameters["CPX_CALLBACK_MIP_DISJCUT"] =      109;
cplex_parameters["CPX_CALLBACK_MIP_FLOWMIR"] =      110;
cplex_parameters["CPX_CALLBACK_MIP_INCUMBENT"] =    111;
cplex_parameters["CPX_CALLBACK_MIP_DELETENODE"] =   112;
cplex_parameters["CPX_CALLBACK_MIP_BRANCH_NOSOLN"] =113;
/* Be sure to check the LP values */

/* Values for getcallbackinfo function */

cplex_parameters["CPX_CALLBACK_INFO_BEST_INTEGER"] =       101;
cplex_parameters["CPX_CALLBACK_INFO_BEST_REMAINING"] =     102;
cplex_parameters["CPX_CALLBACK_INFO_NODE_COUNT"] =         103;
cplex_parameters["CPX_CALLBACK_INFO_NODES_LEFT"] =         104;
cplex_parameters["CPX_CALLBACK_INFO_MIP_ITERATIONS"] =     105;
cplex_parameters["CPX_CALLBACK_INFO_CUTOFF"] =             106;
cplex_parameters["CPX_CALLBACK_INFO_CLIQUE_COUNT"] =       107;
cplex_parameters["CPX_CALLBACK_INFO_COVER_COUNT"] =        108;
cplex_parameters["CPX_CALLBACK_INFO_MIP_FEAS"] =           109;
cplex_parameters["CPX_CALLBACK_INFO_FLOWCOVER_COUNT"] =    110;
cplex_parameters["CPX_CALLBACK_INFO_GUBCOVER_COUNT"] =     111;
cplex_parameters["CPX_CALLBACK_INFO_IMPLBD_COUNT"] =       112;
cplex_parameters["CPX_CALLBACK_INFO_PROBE_PHASE"] =        113;
cplex_parameters["CPX_CALLBACK_INFO_PROBE_PROGRESS"] =     114;
cplex_parameters["CPX_CALLBACK_INFO_FRACCUT_COUNT"] =      115;
cplex_parameters["CPX_CALLBACK_INFO_FRACCUT_PROGRESS"] =   116;
cplex_parameters["CPX_CALLBACK_INFO_DISJCUT_COUNT"] =      117;
cplex_parameters["CPX_CALLBACK_INFO_DISJCUT_PROGRESS"] =   118;
cplex_parameters["CPX_CALLBACK_INFO_FLOWPATH_COUNT"] =     119;
cplex_parameters["CPX_CALLBACK_INFO_MIRCUT_COUNT"] =       120;
cplex_parameters["CPX_CALLBACK_INFO_FLOWMIR_PROGRESS"] =   121;
cplex_parameters["CPX_CALLBACK_INFO_ZEROHALFCUT_COUNT"] =  122;
cplex_parameters["CPX_CALLBACK_INFO_MY_THREAD_NUM"] =      123;
cplex_parameters["CPX_CALLBACK_INFO_USER_THREADS"] =       124;
cplex_parameters["CPX_CALLBACK_INFO_MIP_REL_GAP"] =        125;
cplex_parameters["CPX_CALLBACK_INFO_MCFCUT_COUNT"] =       126;

/* Values for getcallbacknodeinfo function */

cplex_parameters["CPX_CALLBACK_INFO_NODE_SIINF"] =         201;
cplex_parameters["CPX_CALLBACK_INFO_NODE_NIINF"] =         202;
cplex_parameters["CPX_CALLBACK_INFO_NODE_ESTIMATE"] =      203;
cplex_parameters["CPX_CALLBACK_INFO_NODE_DEPTH"] =         204;
cplex_parameters["CPX_CALLBACK_INFO_NODE_OBJVAL"] =        205;
cplex_parameters["CPX_CALLBACK_INFO_NODE_TYPE"] =          206;
cplex_parameters["CPX_CALLBACK_INFO_NODE_VAR"] =           207;
cplex_parameters["CPX_CALLBACK_INFO_NODE_SOS"] =           208;
cplex_parameters["CPX_CALLBACK_INFO_NODE_SEQNUM"] =        209;
cplex_parameters["CPX_CALLBACK_INFO_NODE_USERHANDLE"] =    210;
cplex_parameters["CPX_CALLBACK_INFO_NODE_NODENUM"] =       211;

/* Values for getcallbacksosinfo function */

cplex_parameters["CPX_CALLBACK_INFO_SOS_TYPE"] =           240;
cplex_parameters["CPX_CALLBACK_INFO_SOS_SIZE"] =           241;
cplex_parameters["CPX_CALLBACK_INFO_SOS_IS_FEASIBLE"] =    242;
cplex_parameters["CPX_CALLBACK_INFO_SOS_MEMBER_INDEX"] =   244;
cplex_parameters["CPX_CALLBACK_INFO_SOS_MEMBER_REFVAL"] =  246;
cplex_parameters["CPX_CALLBACK_INFO_SOS_NUM"] =            247;

/* Values for getcallbackindicatorinfo function */

cplex_parameters["CPX_CALLBACK_INFO_IC_NUM"] =             260;
cplex_parameters["CPX_CALLBACK_INFO_IC_IMPLYING_VAR"] =    261;
cplex_parameters["CPX_CALLBACK_INFO_IC_IMPLIED_VAR"] =     262;
cplex_parameters["CPX_CALLBACK_INFO_IC_SENSE"] =           263;
cplex_parameters["CPX_CALLBACK_INFO_IC_COMPL"] =           264;
cplex_parameters["CPX_CALLBACK_INFO_IC_RHS"] =             265;
cplex_parameters["CPX_CALLBACK_INFO_IC_IS_FEASIBLE"] =     266;

/* Value for accessing the incumbent using the solution pool routines */

cplex_parameters["CPX_INCUMBENT_ID"] = -1;

/* Be sure to check the LP values */

/* Callback return codes */

cplex_parameters["CPX_CALLBACK_DEFAULT"] =            0;
cplex_parameters["CPX_CALLBACK_SET"] =                2;
cplex_parameters["CPX_CALLBACK_FAIL"] =               1;

/* For CPXgetnodeintfeas */
cplex_parameters["CPX_INTEGER_FEASIBLE"] =            0;
cplex_parameters["CPX_INTEGER_INFEASIBLE"] =          1;
cplex_parameters["CPX_IMPLIED_INTEGER_FEASIBLE"] =    2;

/* MIP Parameter numbers */
cplex_parameters["CPX_PARAM_BRDIR"] =               2001;
cplex_parameters["CPX_PARAM_BTTOL"] =               2002;
cplex_parameters["CPX_PARAM_CLIQUES"] =             2003;
cplex_parameters["CPX_PARAM_COEREDIND"] =           2004;
cplex_parameters["CPX_PARAM_COVERS"] =              2005;
cplex_parameters["CPX_PARAM_CUTLO"] =               2006;
cplex_parameters["CPX_PARAM_CUTUP"] =               2007;
cplex_parameters["CPX_PARAM_EPAGAP"] =              2008;
cplex_parameters["CPX_PARAM_EPGAP"] =               2009;
cplex_parameters["CPX_PARAM_EPINT"] =               2010;
cplex_parameters["CPX_PARAM_MIPDISPLAY"] =          2012;
cplex_parameters["CPX_PARAM_MIPINTERVAL"] =         2013;
cplex_parameters["CPX_PARAM_INTSOLLIM"] =           2015;
cplex_parameters["CPX_PARAM_NODEFILEIND"] =         2016;
cplex_parameters["CPX_PARAM_NODELIM"] =             2017;
cplex_parameters["CPX_PARAM_NODESEL"] =             2018;
cplex_parameters["CPX_PARAM_OBJDIF"] =              2019;
cplex_parameters["CPX_PARAM_MIPORDIND"] =           2020;
cplex_parameters["CPX_PARAM_RELOBJDIF"] =           2022;
cplex_parameters["CPX_PARAM_STARTALG"] =            2025;
cplex_parameters["CPX_PARAM_SUBALG"] =              2026;
cplex_parameters["CPX_PARAM_TRELIM"] =              2027;
cplex_parameters["CPX_PARAM_VARSEL"] =              2028;
cplex_parameters["CPX_PARAM_BNDSTRENIND"] =         2029;
cplex_parameters["CPX_PARAM_HEURFREQ"] =            2031;
cplex_parameters["CPX_PARAM_MIPORDTYPE"] =          2032;
cplex_parameters["CPX_PARAM_CUTSFACTOR"] =          2033;
cplex_parameters["CPX_PARAM_RELAXPREIND"] =         2034;
cplex_parameters["CPX_PARAM_PRESLVND"] =            2037;
cplex_parameters["CPX_PARAM_BBINTERVAL"] =          2039;
cplex_parameters["CPX_PARAM_FLOWCOVERS"] =          2040;
cplex_parameters["CPX_PARAM_IMPLBD"] =              2041;
cplex_parameters["CPX_PARAM_PROBE"] =               2042;
cplex_parameters["CPX_PARAM_GUBCOVERS"] =           2044;
cplex_parameters["CPX_PARAM_STRONGCANDLIM"] =       2045;
cplex_parameters["CPX_PARAM_STRONGITLIM"] =         2046;
cplex_parameters["CPX_PARAM_FRACCAND"] =            2048;
cplex_parameters["CPX_PARAM_FRACCUTS"] =            2049;
cplex_parameters["CPX_PARAM_FRACPASS"] =            2050;
cplex_parameters["CPX_PARAM_FLOWPATHS"] =           2051;
cplex_parameters["CPX_PARAM_MIRCUTS"] =             2052;
cplex_parameters["CPX_PARAM_DISJCUTS"] =            2053;
cplex_parameters["CPX_PARAM_AGGCUTLIM"] =           2054;
cplex_parameters["CPX_PARAM_MIPCBREDLP"] =          2055    ;
cplex_parameters["CPX_PARAM_CUTPASS"] =             2056;
cplex_parameters["CPX_PARAM_MIPEMPHASIS"] =         2058;
cplex_parameters["CPX_PARAM_SYMMETRY"] =            2059;
cplex_parameters["CPX_PARAM_DIVETYPE"] =            2060;
cplex_parameters["CPX_PARAM_RINSHEUR"] =            2061;
cplex_parameters["CPX_PARAM_SUBMIPNODELIM"] =       2062;
cplex_parameters["CPX_PARAM_LBHEUR"] =              2063;
cplex_parameters["CPX_PARAM_REPEATPRESOLVE"] =      2064;
cplex_parameters["CPX_PARAM_PROBETIME"] =           2065;
cplex_parameters["CPX_PARAM_POLISHTIME"] =          2066;
cplex_parameters["CPX_PARAM_REPAIRTRIES"] =         2067;
cplex_parameters["CPX_PARAM_EPLIN"] =               2068;
cplex_parameters["CPX_PARAM_EPRELAX"] =             2073;
cplex_parameters["CPX_PARAM_FPHEUR"] =              2098;
cplex_parameters["CPX_PARAM_EACHCUTLIM"] =          2102;
cplex_parameters["CPX_PARAM_SOLNPOOLCAPACITY"] =    2103 ;
cplex_parameters["CPX_PARAM_SOLNPOOLREPLACE"] =     2104 ;
cplex_parameters["CPX_PARAM_SOLNPOOLGAP"] =         2105 ;
cplex_parameters["CPX_PARAM_SOLNPOOLAGAP"] =        2106 ;
cplex_parameters["CPX_PARAM_SOLNPOOLINTENSITY"] =   2107;
cplex_parameters["CPX_PARAM_POPULATELIM"] =         2108;
cplex_parameters["CPX_PARAM_MIPSEARCH"] =           2109;
cplex_parameters["CPX_PARAM_MIQCPSTRAT"] =          2110;
cplex_parameters["CPX_PARAM_ZEROHALFCUTS"] =        2111;
cplex_parameters["CPX_PARAM_POLISHAFTEREPAGAP"] =   2126;
cplex_parameters["CPX_PARAM_POLISHAFTEREPGAP"] =    2127;
cplex_parameters["CPX_PARAM_POLISHAFTERNODE"] =     2128;
cplex_parameters["CPX_PARAM_POLISHAFTERINTSOL"] =   2129;
cplex_parameters["CPX_PARAM_POLISHAFTERTIME"] =     2130;
cplex_parameters["CPX_PARAM_MCFCUTS"] =             2134;

/* Solution pools */

/* Values for CPX_PARAM_SOLNPOOLREPLACE */

cplex_parameters["CPX_SOLNPOOL_FIFO"] =   0;

cplex_parameters["CPX_SOLNPOOL_OBJ"] =    1;

cplex_parameters["CPX_SOLNPOOL_DIV"] =    2;


cplex_parameters["CPX_SOLNPOOL_FILTER_DIVERSITY"] =  1;

cplex_parameters["CPX_SOLNPOOL_FILTER_RANGE"] =      2;


/* --------------------------------------------------------------------------
 * File: gcdefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 2005, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

cplex_parameters["CPXERR_UNSUPPORTED_CONSTRAINT_TYPE"] =1212;

cplex_parameters["CPXERR_ILL_DEFINED_PWL"] =1213;


/* Constraint types */

cplex_parameters["CPX_CON_LOWER_BOUND"] =         1;
cplex_parameters["CPX_CON_UPPER_BOUND"] =         2;
cplex_parameters["CPX_CON_LINEAR"] =              3;
cplex_parameters["CPX_CON_QUADRATIC"] =           4;
cplex_parameters["CPX_CON_SOS"] =                 5;
cplex_parameters["CPX_CON_INDICATOR"] =           6;

/* internal types */
cplex_parameters["CPX_CON_MINEXPR"] =             7;
cplex_parameters["CPX_CON_MAXEXPR"] =             8;
cplex_parameters["CPX_CON_PWL"] =                 9;
cplex_parameters["CPX_CON_ABS"] =                 9  /*  same as PWL since using it */;
cplex_parameters["CPX_CON_DISJCST"] =            10;
cplex_parameters["CPX_CON_INDDISJCST"] =         11;
cplex_parameters["CPX_CON_SETVAR"] =             12;
cplex_parameters["CPX_CON_SETVARMEMBER"] =       13;
cplex_parameters["CPX_CON_SETVARCARD"] =         14;
cplex_parameters["CPX_CON_SETVARSUM"] =          15;
cplex_parameters["CPX_CON_SETVARMIN"] =          16;
cplex_parameters["CPX_CON_SETVARMAX"] =          17;
cplex_parameters["CPX_CON_SETVARSUBSET"] =       18;
cplex_parameters["CPX_CON_SETVARDOMAIN"] =       19;
cplex_parameters["CPX_CON_SETVARUNION"] =        20;
cplex_parameters["CPX_CON_SETVARINTERSECTION"] = 21;
cplex_parameters["CPX_CON_SETVARNULLINTERSECT"] =22;
cplex_parameters["CPX_CON_SETVARINTERSECT"] =    23;
cplex_parameters["CPX_CON_SETVAREQ"] =           24;
cplex_parameters["CPX_CON_SETVARNEQ"] =          25;
cplex_parameters["CPX_CON_SETVARNEQCST"] =       26;
cplex_parameters["CPX_CON_LAST_CONTYPE"] =       27;


/* --------------------------------------------------------------------------
 * File: netdefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 1998, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

/* Network parameters */

cplex_parameters["CPX_PARAM_NETITLIM"] =            5001;
cplex_parameters["CPX_PARAM_NETEPOPT"] =            5002;
cplex_parameters["CPX_PARAM_NETEPRHS"] =            5003;
cplex_parameters["CPX_PARAM_NETPPRIIND"] =          5004;
cplex_parameters["CPX_PARAM_NETDISPLAY"] =          5005;


/* NET/MIN files format errors */

cplex_parameters["CPXERR_NET_DATA"] =              1530;
cplex_parameters["CPXERR_NOT_MIN_COST_FLOW"] =     1531;
cplex_parameters["CPXERR_BAD_ROW_ID"] =            1532;
cplex_parameters["CPXERR_BAD_CHAR"] =              1537;
cplex_parameters["CPXERR_NET_FILE_SHORT"] =        1538;


/* NETOPT display values */

cplex_parameters["CPXNET_NO_DISPLAY_OBJECTIVE"] =    0;
cplex_parameters["CPXNET_TRUE_OBJECTIVE"] =          1;
cplex_parameters["CPXNET_PENALIZED_OBJECTIVE"] =     2;

/* NETOPT pricing parameters */

cplex_parameters["CPXNET_PRICE_AUTO"] =              0;
cplex_parameters["CPXNET_PRICE_PARTIAL"] =           1;
cplex_parameters["CPXNET_PRICE_MULT_PART"] =         2;
cplex_parameters["CPXNET_PRICE_SORT_MULT_PART"] =    3;



/* --------------------------------------------------------------------------
 * File: qpdefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 1995, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

/* QP error codes */

cplex_parameters["CPXERR_Q_NOT_POS_DEF"] =  5002;
cplex_parameters["CPXERR_NOT_QP"] =         5004;
cplex_parameters["CPXERR_Q_DUP_ENTRY"] =    5011;
cplex_parameters["CPXERR_Q_NOT_SYMMETRIC"] =5012;
cplex_parameters["CPXERR_Q_NOT_INDEF"] =    5014;

/* Copying data */

cplex_parameters["CPX_PARAM_QPNZREADLIM"] = 4001;

/* Presolve */
cplex_parameters["CPX_PARAM_QPMAKEPSDIND"] =4010;


/* --------------------------------------------------------------------------
 * File: socpdefs.h
 * Version 12.1.0
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5724-Y48
 * (c) Copyright IBM Corporation 2003, 2009. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 *---------------------------------------------------------------------------
 */

/* SOCP error codes */
cplex_parameters["CPXERR_QCP_SENSE"] =      6002;

}
