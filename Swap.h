#pragma once

class Cluster;
class PixelNode;

class Swap
{
public:
	Swap();
	Swap(Cluster* c1, Cluster* c2, PixelNode* px);

	void calcSwapCost();

	// Swaps the px from c1 to c2 changing the energy
	void swap();

public:
	Cluster* c1 = nullptr;
	Cluster* c2 = nullptr;
	PixelNode* px = nullptr;
	float cost = 0.0f;
};