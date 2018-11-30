#include "ClusterPair.h"
#include "Cluster.h"

ClusterPair::ClusterPair(Cluster* c1, Cluster* c2)
{
	ClusterPair::c1 = c1;
	ClusterPair::c2 = c2;
}

// Using the sum of squares shortcut sum(xi^2) + sum(xi)^2/n we can easily avoid resumming
// there is also no need for a mean with the shortcut
void ClusterPair::calcMergingCost()
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

// Returns the node that is not this one
Cluster* ClusterPair::getNeighbor(Cluster* cluster)
{
	if (c1 == cluster)
		return c2;
	else
		return c1;
}