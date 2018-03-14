#pragma once
#include <vector>

// Forward declaration
class ClusterPair;

class PixelNode
{
public:
	PixelNode() { }

	PixelNode(float x, float y, float z, float g)
	{
		PixelNode::x = x;
		PixelNode::y = y;
		PixelNode::z = z;
		PixelNode::g = g;
	}

public:
	float x = -1.0f;
	float y = -1.0f;
	float z = -1.0f;
	float g = -1.0f;
	bool boundaryPx = false;
	Cluster* parent;
};

class Cluster
{
public:
	void calcEnergy()
	{
		unsigned int numPx = static_cast<unsigned int>(pixels.size());
		for (unsigned int i = 0; i < numPx; i++)
		{
			const PixelNode* px = pixels[i];
			sumX += px->x;
			sumY += px->y;
			sumZ += px->z;
			sumG += px->g;
			sumSqr += px->x * px->x + px->y * px->y + px->z * px->z + px->g * px->g;
		}
		energy = sumSqr - (sumX * sumX + sumY * sumY + sumZ * sumZ + sumG * sumG) / numPx;
	}

	void addEdge(ClusterPair* pair)
	{
		// Only add the pair if it isn't already in the pair array
		for (unsigned int i = 0; i < pairs.size(); i++)
		{
			if (pairs[i] == pair)
				return;
		}
		pairs.push_back(pair);
	}

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