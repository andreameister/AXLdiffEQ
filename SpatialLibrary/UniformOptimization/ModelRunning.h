//
//  ModelRunning.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ModelRunning__
#define __UniformOptimization__ModelRunning__

#include "sundials_nvector.h"
#include "cvode_impl.h"
#include "cvode.h"
#include "sundials_dense.h"

#define autocrineT 10000
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

#define Nspecies 13
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define maxR 1.0
#define Nparams 11

struct rates {
    double Binding1;   ///< Forward binding rate for Ig1
    double Binding2;   ///< Forward binding rate for Ig2
    double Unbinding1; ///< Reverse binding rate for Ig1
    double Unbinding2; ///< Reverse binding rate for Ig2
    double xFwd1;      ///< Reaction 1 forward rate.
    double xRev1;      ///< Reaction 1 reverse rate.
    double xFwd2;      ///< Reaction 2 forward rate.
    double xRev2;      ///< Reaction 2 reverse rate.
    double xFwd3;      ///< Reaction 3 forward rate.
    double xRev3;      ///< Reaction 3 reverse rate.
    double xFwd4;      ///< Reaction 4 forward rate.
    double xRev4;      ///< Reaction 4 reverse rate.
    double xFwd5;      ///< Reaction 5 forward rate.
    double xRev5;      ///< Reaction 5 reverse rate.
    double xFwd6;      ///< Reaction 6 forward rate.
    double xRev6;      ///< Reaction 6 reverse rate.
    double expression; ///< AXL expression rate.
    double internalize;///< Non-pY species internalization rate.
    double pYinternalize;///< pY species internalization rate.
    double kRec;       ///< Recycling rate.
    double kDeg;       ///< Degradation rate.
    double fElse;      ///< Recycling fraction for non-D2 species.
    double fD2;        ///< Recycling fraction for D2.
    int pD1;           // Either or 0 or 1, denoting whether the D1 species is active
    double internalFrac;  //Fraction of receptor internalized
    double internalV;  //Internal volume
    double autocrine;  //Autocrine levels of Gas6
    double gasCur;     //Current level of Gas6 (changes based on simulation)
    double *diffD;  
    double *gasProfile;
};

static const double times[2] = {60, 240}; ///< Times of kinetic measurements.
static const double Gass[6] = {64, 16, 4, 1, 0.25, 0}; ///< Kinetic Gas6 doses.
static const double kTPS[5] = {0, 0.5, 1, 5, 10};


static const double pYk[5] = {4.1, 3.6, 5.3, 11.3, 11.6};
static const double pYkErr[5] = {1.0, 1.3, 1.1, 0.7, 0.4};

static const double U87pYk[5] = {0.911, 0.904, 0.935, 1.08, 1.25};
static const double U87pYkErr[5] = {0.134, 0.0991, 0.152, 0.135, 0.130};

static const double U87pY[6][2] = { ///< pY measurements on short time scales.
    {0.00, 0.00},
    {0.89, 1.61},
    {1.74, 1.10},
    {1.35, 1.05},
    {1.93, 1.33},
    {1.00, 1.00}};

static const double U87pYerror[6][2] = { ///< Error for short time scale pY measurements.
    {1000, 1000},
    {0.44, 0.20},
    {0.30, 0.28},
    {0.41, 0.12},
    {0.17, 0.15},
    {0.19, 0.19}};


// Wrapping is outermost cell line, then Gas, then time
static const double pY[6][2] = { ///< pY measurements on short time scales.
    {10.75952427, 8.305264139},
    {7.390399159, 7.056438019},
    {7.144036441, 7.680851079},
    {4.570826833, 8.184089069},
    {6.107714557, 7.204021903},
    {7.535575387, 7.535575387}};

static const double pYerror[6][2] = { ///< Error for short time scale pY measurements.
    {1.74431775,  2.100723242},
    {1.267611,    1.260108508},
    {0.898008437, 1.680415875},
    {1.521844479, 0.871927763},
    {0.932623012, 0.563182873},
    {0.812417951, 0.812417951}};

static const double tot[6][2] = {
    {3443.11, 3219.69},
    {3143.41, 3353.82},
    {3018.88, 3611.82},
    {2608.88, 3448.21},
    {2690.24, 3168.14},
    {2672.00, 2672.00}}; // A549


static const double U87tot[6][2] = {
    {1082.16, 1002},
    {1108.88, 1002},
    {1162.32, 908.48},
    {961.92, 1015.36},
    {1002, 1135.6},
    {1336, 1108.88}};


static const double U87totError[6][2] = {
    {106.88, 53.44},
    {53.44, 40.08},
    {66.8, 146.96},
    {187.04, 120.24},
    {66.8, 66.8},
    {160.32, 66.8}};



static const double totError[6][2] = {
    {174.38, 132.10},
    {189.03, 129.93},
    {245.75, 225.42},
    {154.89, 203.72},
    {128.72, 187.34},
    {82.62, 82.62}}; // A549

static const double surf[6][2] = {
    {0.206, 0.239},
    {0.274, 0.316},
    {0.281, 0.251},
    {0.220, 0.302},
    {0.256, 0.281},
    {0.257, 0.337}}; // A549

static const double surfError[6][2] = {
    {0.043, 0.015},
    {0.047, 0.037},
    {0.032, 0.030},
    {0.025, 0.036},
    {0.044, 0.035},
    {0.030, 0.023}}; // A549


int AXL_react(double, N_Vector, N_Vector, void *); //Takes concentrations, rates of change, and user data to 
double surfCalc (N_Vector);  //Uses to CVode state to calculate amount of receptor on the surface of the cell
double pYcalc (N_Vector, struct rates *);  //Uses CVode state to calculate amount of phosphorylated receptor
double totCalc (const N_Vector, const struct rates * const p);  //Uses CVode state to calculate the total number of receptors
struct rates Param(const double *);  //Creates a construction for reaction constants using a set of input parameters
double calcError (struct rates, double *);  //Calculates the chi squared error for all parts based on comparison of predicted and measured amounts of receptor
double U87calcError (struct rates, double *);  //Performs calcError for the U87 cell line
void*initState(N_Vector, struct rates *);  //Calculates the initial state using CVode going to a steady state
void calcProfile (N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, struct rates *, double, double);  
//Using the outputs of surfCalc, pYcalc, and totCalc and given time points and stimulations uses the CVode solver to create ODEs
void calcProfileSet (double *, double *, double *, double *, double *, struct rates *, unsigned int, double, double *);
void diffusionSolution(double *, double *, unsigned int, double *, double *, unsigned int, double *);
//Models the diffusion of Gas6





#endif /* defined(__UniformOptimization__ModelRunning__) */
