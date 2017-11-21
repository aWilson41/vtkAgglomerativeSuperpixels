#include "vtkSuperpixelFilter.h"

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <algorithm>
#include <time.h>

vtkStandardNewMacro(vtkSuperpixelFilter);

int vtkSuperpixelFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	// Get the info objects
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkImageData* output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	int* dim = input->GetDimensions();
	int numDim = input->GetDataDimension();
	if (numDim != 2 && numDim != 3)
	{
		vtkWarningMacro(<<"Only 2/3d images supported.");
		return 1;
	}
	if (input->GetNumberOfScalarComponents() > 1)
	{
		vtkWarningMacro(<<"Only grayscale images supported.");
		return 1;
	}
	if (numSuperpixels > static_cast<unsigned int>(dim[0] * dim[1] * dim[2]))
	{
		vtkWarningMacro(<<"Number of superpixels larger than number of pixels.");
		return 1;
	}
	if (input->GetScalarType() != 3)
	{
		vtkWarningMacro(<<"Only unsigned char is supported.");
		return 1;
	}
	// Set the information
	output->SetExtent(outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
	output->SetNumberOfScalarComponents(3, outInfo);
	output->SetDimensions(dim);
	output->AllocateScalars(outInfo);
	
	// Get the input/output image pointers
	unsigned char* outPtr = static_cast<unsigned char*>(output->GetScalarPointer());
	unsigned char* inPtr = static_cast<unsigned char*>(input->GetScalarPointer());

	// Start the clock
	std::clock_t start;
	double duration = 0.0;
	start = std::clock();

	// Creates clusters and heap from input image
	if (numDim == 2)
	{
		initClusters(inPtr, dim[0], dim[1], 1);
		initHeap(dim[0], dim[1], 1, numDim);
	}
	else if (numDim == 3)
	{
		initClusters(inPtr, dim[0], dim[1], dim[2]);
		initHeap(dim[0], dim[1], dim[2], numDim);
	}

	unsigned int numPx = dim[0] * dim[1] * dim[2];
	unsigned int n = numPx;
	while (n > numSuperpixels)
	{
		// Pull the pair with the least energy
		ClusterPair* pair = static_cast<ClusterPair*>(minHeap.extract());
		// Merge cluster2 into cluster1
		Cluster* c1 = pair->p1;
		Cluster* c2 = pair->p2;

		// Put all of cluster2's pixels in cluster1
		c1->pixels.insert(c1->pixels.end(), c2->pixels.begin(), c2->pixels.end());
		// The clusters energy, color, and positional info was already calculated during calcMergeCost
		// Now that we are merging we'll just copy that over
		c1->sumX = pair->sumX;
		c1->sumY = pair->sumY;
		c1->sumZ = pair->sumZ;
		c1->sumG = pair->sumG;
		c1->sumSqr = pair->sumSqr;
		c1->energy = pair->energy;

		// Remove bad edges (duplicates/etc)
		removeEdges(pair, c1, c2);

		// Every edge of cluster1 should have energy updated
		for (unsigned int i = 0; i < c1->pairs.size(); i++)
		{
			c1->pairs[i]->calcMergingCost();
			c1->pairs[i]->heap_key(-c1->pairs[i]->dEnergy);
			minHeap.update(c1->pairs[i]);
		}

		// Mark the cluster as removed from the heap (we don't want to have to search for it)
		c2->energy = -1.0;

		// Finally delete the pair
		delete pair;
		n--;
		//printf("%d pixels.\n", static_cast<int>(n));
	}

	// Set random colors to every cluster
	srand(time(NULL));
	for (unsigned int i = 0; i < numPx; i++)
	{
		if (clusters[i].energy != -1.0)
		{
			const unsigned char r = rand() % 255;
			const unsigned char g = rand() % 255;
			const unsigned char b = rand() % 255;
			for (unsigned int j = 0; j < clusters[i].pixels.size(); j++)
			{
				const unsigned int x = clusters[i].pixels[j].x;
				const unsigned int y = clusters[i].pixels[j].y;
				const unsigned int z = clusters[i].pixels[j].z;
				const unsigned int index = (x + (y + dim[1] * z) * dim[0]) * 3;
				outPtr[index] = r;
				outPtr[index + 1] = g;
				outPtr[index + 2] = b;
			}
		}
	}
	delete[] clusters;
	
	// Cleanup the remaining pairs
	for (unsigned int i = 0; i < minHeap.size(); i++)
	{
		delete minHeap.item(i);
	}

	duration = (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	printf("Dim: %d, %d, %d\n", dim[0], dim[1], dim[2]);
	printf("Time: %f\n", duration);
	printf("Superpixels complete\n");

	return 1;
}

void vtkSuperpixelFilter::removeEdges(ClusterPair* pair, Cluster* c1, Cluster* c2)
{
	// Remove current pair from both clusters's pair lists
	c1->pairs.erase(std::remove(c1->pairs.begin(), c1->pairs.end(), pair), c1->pairs.end());
	c2->pairs.erase(std::remove(c2->pairs.begin(), c2->pairs.end(), pair), c2->pairs.end());

	// Duplicate pairs should be removed
	for (unsigned int i = 0; i < c2->pairs.size(); i++)
	{
		if ((c2->pairs[i]->p1 == c1 && c2->pairs[i]->p2 == c2) ||
			(c2->pairs[i]->p1 == c2 && c2->pairs[i]->p2 == c1))
		{
			// Remove it from the minheap (this requires a search of O(n) that could be eliminated through further programming)
			minHeap.remove(c2->pairs[i]);
			// Remove the pair from both node's edge lists
			c2->pairs.erase(c2->pairs.begin() + i);
			i--;
		}
	}
	// Duplicate pairs should be removed
	for (unsigned int i = 0; i < c1->pairs.size(); i++)
	{
		if ((c1->pairs[i]->p1 == c1 && c1->pairs[i]->p2 == c2) ||
			(c1->pairs[i]->p1 == c2 && c1->pairs[i]->p2 == c1))
		{
			// Remove it from the minheap (this requires a search of O(n) that could be eliminated through further programming)
			minHeap.remove(c1->pairs[i]);
			// Remove the pair from both node's edge lists
			c1->pairs.erase(c1->pairs.begin() + i);
			i--;
		}
	}
	// Every pair in cluster2 that is left is redirect to cluster1
	for (unsigned int i = 0; i < c2->pairs.size(); i++)
	{
		if (c2->pairs[i]->p1 == c2)
			c2->pairs[i]->p1 = c1;
		else if (c2->pairs[i]->p2 == c2)
			c2->pairs[i]->p2 = c1;
	}

	// Merge the remaining pairs from cluster2 to cluster1
	c1->pairs.insert(c1->pairs.end(), c2->pairs.begin(), c2->pairs.end());
}

void vtkSuperpixelFilter::initClusters(unsigned char* inPtr, int width, int height, int depth)
{
	//clusters = std::vector<Cluster*>(width * height * depth);
	clusters = new Cluster[width * height * depth];
	int index = 0;
	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				// Create a new node
				PixelNode node(x, y, z, inPtr[index]);
				// Each cluster is a collection of PixelNodes we start with one cluster for every PixelNode then merge clusters
				clusters[index] = Cluster();
				// Add the pixel to the cluster
				clusters[index].pixels.push_back(node);
				clusters[index].calcEnergy();
				index++;
			}
		}
	}
}

void vtkSuperpixelFilter::initHeap(int width, int height, int depth, int numDim)
{
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
				minHeap.insert(pair, -pair->dEnergy);
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
				minHeap.insert(pair, -pair->dEnergy);
			}
		}
	}

	// Depth is 1 so for 2d images we don't want to duplicate edges
	if (numDim == 3)
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
					minHeap.insert(pair, -pair->dEnergy);
				}
			}
		}
	}
}