#pragma once
#include "Cluster.h"
#include "Mx\MxHeap.h"

class ClusterPair : public MxHeapable
{
public:
	ClusterPair() { }

	ClusterPair(Cluster* c1, Cluster* c2)
	{
		ClusterPair::c1 = c1;
		ClusterPair::c2 = c2;
	}

	// Using the sum of squares shortcut sum(xi^2) + sum(xi)^2/n we can easily avoid resumming
	// there is also no need for a mean with the shortcut
	void calcMergingCost()
	{
		// Calculate sums
		sumX = c1->sumX + c2->sumX;
		sumY = c1->sumY + c2->sumY;
		sumZ = c1->sumZ + c2->sumZ;
		sumG = c1->sumG + c2->sumG;
		sumSqr = c1->sumSqr + c2->sumSqr;

		// Calculate the resulting energy and change in energy
		energy = sumSqr - (sumX * sumX + sumY * sumY + sumZ * sumZ + sumG * sumG) / (c1->pixels.size() + c2->pixels.size());
		dEnergy = energy - c1->energy - c2->energy;
	}

public:
	Cluster* c1;
	Cluster* c2;
	// Merging cost
	float dEnergy = 0.0f;
	float energy = 0.0f;
	float sumX = 0.0f;
	float sumY = 0.0f;
	float sumZ = 0.0f;
	float sumG = 0.0f;
	float sumSqr = 0.0f;
};