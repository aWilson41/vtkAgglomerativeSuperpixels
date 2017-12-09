#pragma once

#include "Cluster.h"
#include "Mx\MxHeap.h"

class ClusterPair : public MxHeapable
{
public:
	ClusterPair() { }

	ClusterPair(Cluster* p1, Cluster* p2)
	{
		ClusterPair::p1 = p1;
		ClusterPair::p2 = p2;
	}

	// Using the sum of squares shortcut sum(xi^2) + sum(xi)^2/n we can easily avoid resumming
	// there is also no need for a mean with the shortcut
	void calcMergingCost()
	{
		// Calculate sums
		sumX = p1->sumX + p2->sumX;
		sumY = p1->sumY + p2->sumY;
		sumZ = p1->sumZ + p2->sumZ;
		sumG = p1->sumG + p2->sumG;
		sumSqr = p1->sumSqr + p2->sumSqr;

		// Calculate the resulting energy and change in energy
		energy = sumSqr - (sumX * sumX + sumY * sumY + sumZ * sumZ + sumG * sumG) / (p1->pixels.size() + p2->pixels.size());
		dEnergy = energy - p1->energy - p2->energy;
	}

public:
	Cluster* p1;
	Cluster* p2;
	// Merging cost
	float dEnergy = 0.0;
	float energy = 0.0;
	float sumX = 0.0;
	float sumY = 0.0;
	float sumZ = 0.0;
	float sumG = 0.0;
	float sumSqr = 0.0;
};