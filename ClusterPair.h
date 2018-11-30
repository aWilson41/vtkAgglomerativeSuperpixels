#pragma once
#include "Mx/MxHeap.h"

class Cluster;

class ClusterPair : public MxHeapable
{
public:
	ClusterPair() { }

	ClusterPair(Cluster* c1, Cluster* c2);

	void calcMergingCost();

	Cluster* getNeighbor(Cluster* cluster);

public:
	Cluster* c1 = nullptr;
	Cluster* c2 = nullptr;
	// Merging cost
	float dEnergy = 0.0f;
	float energy = 0.0f;
	float sumX = 0.0f;
	float sumY = 0.0f;
	float sumZ = 0.0f;
	float sumG = 0.0f;
	float sumSqr = 0.0f;
};
