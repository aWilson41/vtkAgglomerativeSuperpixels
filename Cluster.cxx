#include "Cluster.h"
#include "ClusterPair.h"
#include "PixelNode.h"

void Cluster::calcEnergy()
{
	unsigned int numPx = static_cast<unsigned int>(pixels.size());
	sumX = sumY = sumZ = sumG = sumSqr = 0.0f;
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

void Cluster::addEdge(ClusterPair* pair)
{
	// Only add the pair if it isn't already in the pair array
	for (unsigned int i = 0; i < pairs.size(); i++)
	{
		if (pairs[i] == pair)
			return;
	}
	pairs.push_back(pair);
}

float* Cluster::getCentroid()
{
	static float c[3] = {
		sumX / static_cast<float>(pixels.size()),
		sumY / static_cast<float>(pixels.size()),
		sumZ / static_cast<float>(pixels.size()) };
	return c;
}
float Cluster::getAvgIntensity() { return sumG / static_cast<float>(pixels.size()); }
float Cluster::getMaxIntensity()
{
	float max = std::numeric_limits<float>::min();
	for (unsigned int i = 0; i < pixels.size(); i++)
	{
		if (pixels[i]->g > max)
			max = pixels[i]->g;
	}
	return max;
}
float Cluster::getMinIntensity()
{
	float min = std::numeric_limits<float>::max();
	for (unsigned int i = 0; i < pixels.size(); i++)
	{
		if (pixels[i]->g < min)
			min = pixels[i]->g;
	}
	return min;
}