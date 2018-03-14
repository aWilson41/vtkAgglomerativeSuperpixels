#include "vtkSuperpixelFilter.h"
#include "ClusterPair.h"
#include "Cluster.h"
#include "Swap.h"

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <algorithm>

// For benchmarking
#include <chrono>

vtkStandardNewMacro(vtkSuperpixelFilter);

int vtkSuperpixelFilter::RequestInformation(vtkInformation* vtkNotUsed(request), vtkInformationVector** vtkNotUsed(inputVec), vtkInformationVector* outputVec)
{
	// get the info objects
	vtkInformation* outInfo = outputVec->GetInformationObject(0);
	if (outputType == AVGCOLOR || outputType == LABEL)
		vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
	else if (outputType == RANDRGB)
		vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_UNSIGNED_CHAR, 3);
	return 1;
}

int vtkSuperpixelFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVec, vtkInformationVector* outputVec)
{
	// Get input image
	vtkInformation* inInfo = inputVec[0]->GetInformationObject(0);
	vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	// Get output
	vtkInformation* outInfo = outputVec->GetInformationObject(0);
	vtkImageData* output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	int* dim = input->GetDimensions();
	int numDim = input->GetDataDimension();
	if (numDim != 2 && numDim != 3)
	{
		vtkWarningMacro(<< "Only 2/3d images supported.");
		return 1;
	}
	if (input->GetNumberOfScalarComponents() > 1)
	{
		vtkWarningMacro(<< "Only grayscale images supported.");
		return 1;
	}
	unsigned int numPx = dim[0] * dim[1] * dim[2];
	if (numSuperpixels >= numPx || numSuperpixels <= 0)
	{
		vtkWarningMacro(<< "Invalid number of superpixles.");
		return 1;
	}
	if (input->GetScalarType() != VTK_FLOAT)
	{
		vtkWarningMacro(<< "Only float images are supported.");
		return 1;
	}
	// Set the information
	output->SetExtent(outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
	output->SetNumberOfScalarComponents(1, outInfo);
	output->SetScalarType(VTK_FLOAT, outInfo);
	if (outputType == RANDRGB)
	{
		output->SetNumberOfScalarComponents(3, outInfo);
		output->SetScalarType(VTK_UNSIGNED_CHAR, outInfo);
	}
	output->SetDimensions(dim);
	output->AllocateScalars(outInfo);

	// Creates clusters and heap from input image
	initClusters(input);
	initHeap(input);

	auto start = std::chrono::steady_clock::now();

	unsigned int n = numPx;
	// For progress
	unsigned int pxToDecimate = numPx - numSuperpixels;
	// If exclude 0 is on it may try to continue without anything in the heap
	while (n > numSuperpixels)
	{
		// Pull the pair with the least energy
		ClusterPair* pair = static_cast<ClusterPair*>(minHeap->extract());
		// Merge cluster2 into cluster1
		Cluster* c1 = pair->c1;
		Cluster* c2 = pair->c2;

		// Put all of cluster2's pixels in cluster1
		c1->pixels.insert(c1->pixels.end(), c2->pixels.begin(), c2->pixels.end());
		// The clusters energy, color, and positional info was already calculated before putting in minheap
		// Now that we are merging we'll just copy that over
		c1->sumX = pair->sumX;
		c1->sumY = pair->sumY;
		c1->sumZ = pair->sumZ;
		c1->sumG = pair->sumG;
		c1->sumSqr = pair->sumSqr;
		c1->energy = pair->energy;

		// Remove bad edges (duplicates/etc)
		removeEdges(pair);

		// Every edge of cluster1 should have energy updated
		for (unsigned int i = 0; i < c1->pairs.size(); i++)
		{
			c1->pairs[i]->calcMergingCost();
			c1->pairs[i]->heap_key(-c1->pairs[i]->dEnergy);
			minHeap->update(c1->pairs[i]);
		}

		// Mark the cluster as removed from the heap
		c2->energy = -1.0f;

		// Finally delete the pair
		delete pair;
		n--;

		// Update the progress of the filter
		UpdateProgress((numPx - n) / static_cast<double>(pxToDecimate));
	}

	// Extract a vector of clusters that are still valid
	finalClusters.clear();
	finalClusters.resize(numSuperpixels);
	unsigned int index = 0;
	for (unsigned int i = 0; i < numPx; i++)
	{
		if (clusters[i].energy != -1.0f)
			finalClusters[index++] = &clusters[i];
	}

	if (swap)
	{
		for (unsigned int i = 0; i < swapIterations; i++)
		{
			computeSwap(dim[0], dim[1], dim[2]);
		}
	}

	// User can specify different output options
	if (outputType == AVGCOLOR)
		calcAvgColors(output);
	else if (outputType == LABEL)
		calcColorLabels(output);
	else if (outputType == RANDRGB)
		calcRandRgb(output);

	// Cleanup
	delete[] clusters;
	delete[] px;
	for (unsigned int i = 0; i < minHeap->size(); i++)
	{
		delete minHeap->item(i);
	}
	delete minHeap;

	auto end = std::chrono::steady_clock::now();
	printf("Time: %f\n", std::chrono::duration<double, std::milli>(end - start).count() / 1000.0);

	return 1;
}

void vtkSuperpixelFilter::removeEdges(ClusterPair* pair)
{
	// We will take all of c2's pairs and add them to c1. But first we need to remove a few pairs
	Cluster* c1 = pair->c1;
	Cluster* c2 = pair->c2;

	// Remove the current pair from both clusters pair list
	c1->pairs.erase(std::remove(c1->pairs.begin(), c1->pairs.end(), pair), c1->pairs.end());
	c2->pairs.erase(std::remove(c2->pairs.begin(), c2->pairs.end(), pair), c2->pairs.end());

	// Redirect every pair from cluster2 to cluster1
	for (int i = 0; i < c2->pairs.size(); i++)
	{
		ClusterPair* pair1 = c2->pairs[i];
		// We know c2 is one of the clusters in the pair. We don't know which one.
		// If p1 is c2, then change the connection so that it's between c1
		if (pair1->c1 == c2)
			pair1->c1 = c1;
		// If not then p2 is c2, change the connection so that it's between c1
		else
			pair1->c2 = c1;
	}

	// Insert all of cluster2's pairs into cluster1
	c1->pairs.insert(c1->pairs.end(), c2->pairs.begin(), c2->pairs.end());

	// Compare every pair in c1 with every other pair in the cluster to remove duplicates
	for (int i = 0; i < c1->pairs.size(); i++)
	{
		ClusterPair* pair1 = c1->pairs[i];
		for (int j = i + 1; j < c1->pairs.size(); j++)
		{
			ClusterPair* pair2 = c1->pairs[j];
			// Remove if they point to the same two clusters
			if ((pair1->c1 == pair2->c1 && pair1->c2 == pair2->c2) ||
				(pair1->c1 == pair2->c2 && pair1->c2 == pair2->c1))
			{
				minHeap->remove(pair2);

				// Not only do we need to remove the pair from the minHeap we need to remove the pair from both of it's clusters pair lists
				c1->pairs.erase(c1->pairs.begin() + j);
				// If p1 is cl, then we should remove the pair from p2's pair list
				if (pair2->c1 == c1)
					pair2->c2->pairs.erase(std::remove(pair2->c2->pairs.begin(), pair2->c2->pairs.end(), pair2), pair2->c2->pairs.end());
				else // If not, then p2 is c1, so we remove the pair from p1's pair list
					pair2->c1->pairs.erase(std::remove(pair2->c1->pairs.begin(), pair2->c1->pairs.end(), pair2), pair2->c1->pairs.end());
				j--;
			}
		}
	}
}

void vtkSuperpixelFilter::computeSwap(int width, int height, int depth)
{
	// Setup parentage for the clusters so every pixel in the image knows what cluster it is apart of
	for (unsigned int i = 0; i < finalClusters.size(); i++)
	{
		for (unsigned int j = 0; j < finalClusters[i]->pixels.size(); j++)
		{
			finalClusters[i]->pixels[j]->parent = finalClusters[i];
		}
		finalClusters[i]->calcEnergy();
	}

	std::vector<Swap> swaps;

	// Neighborhood iteration
	const int neighborhood[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
	//int a = width * height;
	const int neighborShift[4] = { 1, -1, width, -width };
	// Calc the cost of swapping every border pixel to it's axial neighbors
	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				// Find a neighbor cluster to swap this pixel too with the min cost
				int i = x + (y + height * z) * width;
				PixelNode* swapPx = &px[i];
				Cluster* c1 = swapPx->parent;
				Swap minSwap;
				// For every neighbor
				for (unsigned int j = 0; j < 4; j++)
				{
					int xp = x + neighborhood[j][0];
					int yp = y + neighborhood[j][1];
					// Check bounds
					if (xp >= 0 && yp >= 0 && xp < width && yp < height)
					{
						PixelNode* neighborPx = &px[i + neighborShift[j]];
						Cluster* c2 = neighborPx->parent;
						// If the two pixels belong too different clusters
						if (c1 != c2)
						{
							// Create a new potential swap
							Swap newSwap = Swap(c1, c2, swapPx);
							newSwap.calcSwapCost();
							// If the cost is smaller and the px we want to swap hasn't already been reassigned
							if (newSwap.cost < minSwap.cost)
								minSwap = newSwap;
						}
					}
				}
				// As long as the min swap reduces the energy then it's valid
				if (minSwap.cost < 0.0f)
					swaps.push_back(minSwap);
			}
		}
	}
	// Do all the swaps
	for (int i = 0; i < swaps.size(); i++)
	{
		swaps[i].swap();
	}
}


void vtkSuperpixelFilter::calcColorLabels(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster
	float g = 0.0f;
	for (int i = 0; i < finalClusters.size(); i++)
	{
		Cluster* cluster = finalClusters[i];
		// Color every pixel in this cluster this label
		for (unsigned int j = 0; j < cluster->pixels.size(); j++)
		{
			const unsigned int x = cluster->pixels[j]->x;
			const unsigned int y = cluster->pixels[j]->y;
			const unsigned int z = cluster->pixels[j]->z;
			const unsigned int index = x + (y + dim[1] * z) * dim[0];
			outPtr[index] = g;
		}
		// Increment the label
		g += 1.0f;
	}
}

void vtkSuperpixelFilter::calcAvgColors(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster set a unique color
	for (int i = 0; i < finalClusters.size(); i++)
	{
		Cluster* cluster = finalClusters[i];
		float avg = cluster->sumG / cluster->pixels.size();
		for (unsigned int j = 0; j < cluster->pixels.size(); j++)
		{
			const unsigned int x = cluster->pixels[j]->x;
			const unsigned int y = cluster->pixels[j]->y;
			const unsigned int z = cluster->pixels[j]->z;
			const unsigned int index = x + (y + dim[1] * z) * dim[0];
			outPtr[index] = avg;
		}
	}
}

void vtkSuperpixelFilter::calcRandRgb(vtkImageData* output)
{
	unsigned char* outPtr = static_cast<unsigned char*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster set a random rgb color
	for (int i = 0; i < finalClusters.size(); i++)
	{
		Cluster* cluster = finalClusters[i];
		// Create a random color and assign it to every pixel in the cluster
		float r = static_cast<unsigned char>(rand() % 255);
		float g = static_cast<unsigned char>(rand() % 255);
		float b = static_cast<unsigned char>(rand() % 255);
		for (unsigned int j = 0; j < cluster->pixels.size(); j++)
		{
			const unsigned int x = cluster->pixels[j]->x;
			const unsigned int y = cluster->pixels[j]->y;
			const unsigned int z = cluster->pixels[j]->z;
			const unsigned int index = (x + (y + dim[1] * z) * dim[0]) * 3;
			outPtr[index] = r;
			outPtr[index + 1] = g;
			outPtr[index + 2] = b;
		}
	}
}


// Sets up all the clusters and calculates their energy. Returns the number of zeros found in the image.
void vtkSuperpixelFilter::initClusters(vtkImageData* input)
{
	float* inPtr = static_cast<float*>(input->GetScalarPointer());
	int* dim = input->GetDimensions();
	int numPx = dim[0] * dim[1] * dim[2];
	clusters = new Cluster[numPx];
	px = new PixelNode[numPx];
	int index = 0;
	for (int z = 0; z < dim[2]; z++)
	{
		for (int y = 0; y < dim[1]; y++)
		{
			for (int x = 0; x < dim[0]; x++)
			{
				// Each cluster is a collection of PixelNodes we start with one cluster for every PixelNode then merge clusters
				px[index] = PixelNode(x, y, z, inPtr[index] * colorWeight);
				clusters[index] = Cluster();
				// Add the pixel to the cluster
				clusters[index].pixels.push_back(&px[index]);
				clusters[index].calcEnergy();
				index++;
			}
		}
	}
}

// Adds the pairs to the min heap
void vtkSuperpixelFilter::initHeap(vtkImageData* input)
{
	minHeap = new MxHeap();
	int* dim = input->GetDimensions();
	int width = dim[0];
	int height = dim[1];
	int depth = dim[2];
	// Add all the horz edges
	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width - 1; x++)
			{
				int index = x + (y + height * z) * width;
				int index1 = index + 1;
				ClusterPair* pair = new ClusterPair(&clusters[index], &clusters[index1]);
				clusters[index].addEdge(pair);
				clusters[index1].addEdge(pair);
				pair->calcMergingCost();
				minHeap->insert(pair, -pair->dEnergy);
			}
		}
	}
	// Add all the vert edges
	for (int z = 0; z < depth; z++)
	{
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height - 1; y++)
			{
				int index = x + (y + height * z) * width;
				int index1 = index + width;
				ClusterPair* pair = new ClusterPair(&clusters[index], &clusters[index1]);
				clusters[index].addEdge(pair);
				clusters[index1].addEdge(pair);
				pair->calcMergingCost();
				minHeap->insert(pair, -pair->dEnergy);
			}
		}
	}

	// If depth is 1 then we are working with a 2d image and shouldn't add these pairs
	if (depth != 1)
	{
		// Add all the vert edges
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				for (int z = 0; z < depth - 1; z++)
				{
					int index = x + (y + height * z) * width;
					int index1 = index + width * height;
					ClusterPair* pair = new ClusterPair(&clusters[index], &clusters[index1]);
					clusters[index].addEdge(pair);
					clusters[index1].addEdge(pair);
					pair->calcMergingCost();
					minHeap->insert(pair, -pair->dEnergy);
				}
			}
		}
	}
}