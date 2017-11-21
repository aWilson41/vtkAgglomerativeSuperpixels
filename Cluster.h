#pragma once

#include <vector>

// Forward declaration
class ClusterPair;

class PixelNode
{
public:
	PixelNode()
	{
		x = y = z = g = -1.0;
	}

	PixelNode(double x, double y, double z, double g)
	{
		PixelNode::x = x;
		PixelNode::y = y;
		PixelNode::z = z;
		PixelNode::g = g;
	}

public:
	double x, y, z, g;
	double x2, y2, z2, g2;
};

class Cluster
{
public:
	void calcEnergy()
	{
		unsigned int numPx = static_cast<unsigned int>(pixels.size());
		double invSize = 1.0 / numPx;
		for (unsigned int i = 0; i < numPx; i++)
		{
			const PixelNode* px = &pixels[i];
			sumX += px->x;
			sumY += px->y;
			sumZ += px->z;
			sumG += px->g;
			sumSqr += px->x * px->x + px->y * px->y + px->z * px->z + px->g * px->g;
		}
		energy = sumSqr - (sumX * sumX + sumY * sumY + sumZ * sumZ + sumG * sumG) * invSize;
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
	std::vector<PixelNode> pixels;
	// All the edges connected to this cluster
	std::vector<ClusterPair*> pairs;
	// double sum of x, y, z, g components and sqr sum
	double sumX = 0.0;
	double sumY = 0.0;
	double sumZ = 0.0;
	double sumG = 0.0;
	double sumSqr = 0.0;
	// Energy (Spatial + Color)
	double energy = 0.0;
};