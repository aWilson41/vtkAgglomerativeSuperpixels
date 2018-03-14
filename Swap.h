#pragma once
#include "Cluster.h"
#include <algorithm>

class Swap
{
public:
	Swap() { }
	Swap(Cluster* c1, Cluster* c2, PixelNode* px)
	{
		Swap::c1 = c1;
		Swap::c2 = c2;
		Swap::px = px;
	}

	void calcSwapCost()
	{
		float pxSumSqr = px->x * px->x + px->y * px->y + px->z * px->z + px->g * px->g;
		// The energy of c1 if we removed px
		float postSumXC1 = c1->sumX - px->x;
		float postSumYC1 = c1->sumY - px->y;
		float postSumZC1 = c1->sumZ - px->z;
		float postSumGC1 = c1->sumG - px->g;
		float postSumSqrC1 = c1->sumSqr - pxSumSqr;
		float postEnergyC1 = postSumSqrC1 - (postSumXC1 * postSumXC1 + postSumYC1 * postSumYC1 + postSumZC1 * postSumZC1 + postSumGC1 * postSumGC1) / (c1->pixels.size() - 1);
		// The energy of c2 if we added px
		float postSumXC2 = c2->sumX + px->x;
		float postSumYC2 = c2->sumY + px->y;
		float postSumZC2 = c2->sumZ + px->z;
		float postSumGC2 = c2->sumG + px->g;
		float postSumSqrC2 = c2->sumSqr + pxSumSqr;
		float postEnergyC2 = postSumSqrC2 - (postSumXC2 * postSumXC2 + postSumYC2 * postSumYC2 + postSumZC2 * postSumZC2 + postSumGC2 * postSumGC2) / (c2->pixels.size() + 1);

		cost = postEnergyC1 + postEnergyC2 - c1->energy - c2->energy;
	}

	// Swaps the px from c1 to c2 changing the energy
	void swap()
	{
		// Remove pointer to swapPx from c1
		c1->pixels.erase(std::remove(c1->pixels.begin(), c1->pixels.end(), px), c1->pixels.end());
		// Add swapPx to c2
		c2->pixels.push_back(px);
	}

public:
	Cluster * c1;
	Cluster* c2;
	PixelNode* px;
	float cost = std::numeric_limits<float>::max();
};