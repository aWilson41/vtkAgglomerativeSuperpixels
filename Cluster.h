#pragma once
#include <vector>

// Forward declaration
class ClusterPair;
class PixelNode;

class Cluster
{
public:
	void calcEnergy();

	void addEdge(ClusterPair* pair);

	float* getCentroid();
	float getAvgIntensity();
	float getMaxIntensity();
	float getMinIntensity();

public:
	// All the pixels in this cluster
	std::vector<PixelNode*> pixels;
	// All the edges connected to this cluster
	std::vector<ClusterPair*> pairs;
	// double sum of x, y, z, g components and sqr sum
	float sumX = 0.0f;
	float sumY = 0.0f;
	float sumZ = 0.0f;
	float sumG = 0.0f;
	float sumSqr = 0.0f;
	// Energy (Spatial + Color)
	float energy = 0.0f;
};