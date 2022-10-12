#ifndef TESTBENCH_H
#define TESTBENCH_H

// Definitions
#include "myLib.h"
#include "Utilities.h"
#include "Bessel.h"
#include "QuadL.h"
#include "Engine.h"
#include "MLayers.h"

// Functions
void TestPaulus();
void TestChew1();
void TestChew2();
void TestFEKO1();
void TestFEKO2();
void TestReflection(Config *myConfig);
void TestTLGF(Config *myConfig, double z_);
void TestTLGFr(Config *myConfig, double z_);
void TestGoldKretschmann();

#endif