#include "vtkSuperpixelFilter.h"

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <algorithm>

vtkStandardNewMacro(vtkSuperpixelFilter);

int vtkSuperpixelFilter::RequestInformation(vtkInformation* vtkNotUsed(request), vtkInformationVector** vtkNotUsed(inputVec), vtkInformationVector* outputVec)
{
	// get the info objects
	vtkInformation* outInfo = outputVec->GetInformationObject(0);
	if (outputType == AVGCOLOR || outputType == LABEL)
		vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
	else if (outputType == RANDRGB)
		vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 3);
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
	if (numSuperpixels > numPx)
	{
		vtkWarningMacro(<< "Number of superpixels larger than number of pixels.");
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
	if (outputType == RANDRGB)
		output->SetNumberOfScalarComponents(3, outInfo);
	output->SetScalarType(VTK_FLOAT, outInfo);
	output->SetDimensions(dim);
	output->AllocateScalars(outInfo);

	// Creates clusters and heap from input image
	unsigned int zeroCount = initClusters(input);
	if (excludeZero)
	{
		initHeapExcludeZero(input);
		numPx -= zeroCount; // This is how many clusters there are, not including zeros
	}
	else
		initHeap(input);

	std::clock_t start;
	double duration = 0.0;
	start = std::clock();

	unsigned int n = numPx;
	// For progress
	unsigned int pxToDecimate = numPx - numSuperpixels;
	// If exclude 0 is on it may try to continue without anything in the heap
	while (n > numSuperpixels && minHeap.size() > 0)
	{
		// Pull the pair with the least energy
		ClusterPair* pair = static_cast<ClusterPair*>(minHeap.extract());
		// Merge cluster2 into cluster1
		Cluster* c1 = pair->p1;
		Cluster* c2 = pair->p2;

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
			minHeap.update(c1->pairs[i]);
		}

		// Mark the cluster as removed from the heap
		c2->energy = -1.0f;

		// Finally delete the pair
		delete pair;
		n--;

		// Update the progress of the filter
		UpdateProgress((numPx - n) / static_cast<double>(pxToDecimate));
	}

	// User has option to output average color of clusters or labeled clusters
	if (outputType == AVGCOLOR)
		calcAvgColors(output);
	else if (outputType == LABEL)
	{
		if (excludeZero)
			calcColorLabelsExcludeZero(output);
		else
			calcColorLabels(output);
	}
	else if (outputType == RANDRGB)
	{
		if (excludeZero)
			calcRandRgbExcludeZero(output);
		else
			calcRandRgb(output);
	}

	// Cleanup
	delete[] clusters;
	for (unsigned int i = 0; i < minHeap.size(); i++)
	{
		delete minHeap.item(i);
	}

	duration = (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	printf("Time: %f\n", duration);

	return 1;
}

void vtkSuperpixelFilter::removeEdges(ClusterPair* pair)
{
	// We will take all of c2's pairs and add them to c1. But first we need to remove a few pairs
	Cluster* c1 = pair->p1;
	Cluster* c2 = pair->p2;

	// Remove the current pair from both clusters pair list
	c1->pairs.erase(std::remove(c1->pairs.begin(), c1->pairs.end(), pair), c1->pairs.end());
	c2->pairs.erase(std::remove(c2->pairs.begin(), c2->pairs.end(), pair), c2->pairs.end());

	// Redirect every pair from cluster2 to cluster1
	for (int i = 0; i < c2->pairs.size(); i++)
	{
		ClusterPair* pair1 = c2->pairs[i];
		// We know c2 is one of the clusters in the pair. We don't know which one.
		// If p1 is c2, then change the connection so that it's between c1
		if (pair1->p1 == c2)
			pair1->p1 = c1;
		// If not then p2 is c2, change the connection so that it's between c1
		else
			pair1->p2 = c1;
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
			if ((pair1->p1 == pair2->p1 && pair1->p2 == pair2->p2) ||
				(pair1->p1 == pair2->p2 && pair1->p2 == pair2->p1))
			{
				minHeap.remove(pair2);

				// Not only do we need to remove the pair from the minHeap we need to remove the pair from both of it's clusters pair lists
				c1->pairs.erase(c1->pairs.begin() + j);
				// If p1 is cl, then we should remove the pair from p2's pair list
				if (pair2->p1 == c1)
					pair2->p2->pairs.erase(std::remove(pair2->p2->pairs.begin(), pair2->p2->pairs.end(), pair2), pair2->p2->pairs.end());
				else // If not, then p2 is c1, so we remove the pair from p1's pair list
					pair2->p1->pairs.erase(std::remove(pair2->p1->pairs.begin(), pair2->p1->pairs.end(), pair2), pair2->p1->pairs.end());
				j--;
			}
		}
	}
}


void vtkSuperpixelFilter::calcColorLabels(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster
	float g = 0.0f;
	for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++)
	{
		Cluster* cluster = &clusters[i];
		// As long as its energy isn't -1 it is still in the heap
		if (cluster->energy != -1.0f)
		{
			// Color every pixel in this cluster this label
			for (unsigned int j = 0; j < cluster->pixels.size(); j++)
			{
				const unsigned int x = cluster->pixels[j].x;
				const unsigned int y = cluster->pixels[j].y;
				const unsigned int z = cluster->pixels[j].z;
				const unsigned int index = x + (y + dim[1] * z) * dim[0];
				outPtr[index] = g;
			}
			// Increment the label
			g += 1.0f;
		}
	}
}

void vtkSuperpixelFilter::calcColorLabelsExcludeZero(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	float g = 1.0f;
	// For every cluster
	for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++)
	{
		Cluster* cluster = &clusters[i];
		// As long as its energy isn't -1 it is still in the heap
		if (cluster->energy != -1.0f)
		{
			// If G is over 0 label it
			if (cluster->sumG > 0.0f)
			{
				// Color every pixel in this cluster this label
				for (unsigned int j = 0; j < cluster->pixels.size(); j++)
				{
					const unsigned int x = cluster->pixels[j].x;
					const unsigned int y = cluster->pixels[j].y;
					const unsigned int z = cluster->pixels[j].z;
					const unsigned int index = x + (y + dim[1] * z) * dim[0];
					outPtr[index] = g;
				}
				// Increment the label
				g += 1.0f;
			}
			else
			{
				// If not just set the color of the cluster to 0 since it was never used
				const unsigned int x = cluster->pixels[0].x;
				const unsigned int y = cluster->pixels[0].y;
				const unsigned int z = cluster->pixels[0].z;
				const unsigned int index = x + (y + dim[1] * z) * dim[0];
				outPtr[index] = 0.0f;
			}
		}
	}
}

void vtkSuperpixelFilter::calcAvgColors(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster set a unique color
	for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++)
	{
		Cluster* cluster = &clusters[i];
		if (cluster->energy != -1.0f)
		{
			float avg = cluster->sumG / cluster->pixels.size();
			for (unsigned int j = 0; j < cluster->pixels.size(); j++)
			{
				const unsigned int x = cluster->pixels[j].x;
				const unsigned int y = cluster->pixels[j].y;
				const unsigned int z = cluster->pixels[j].z;
				const unsigned int index = x + (y + dim[1] * z) * dim[0];
				outPtr[index] = avg;
			}
		}
	}
}

void vtkSuperpixelFilter::calcRandRgb(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster set a random rgb color
	for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++)
	{
		Cluster* cluster = &clusters[i];
		// If it is still in the heap
		if (cluster->energy != -1.0f)
		{
			// Create a random color and assign it to every pixel in the cluster
			float r = static_cast<float>(rand() % 255);
			float g = static_cast<float>(rand() % 255);
			float b = static_cast<float>(rand() % 255);
			for (unsigned int j = 0; j < cluster->pixels.size(); j++)
			{
				const unsigned int x = cluster->pixels[j].x;
				const unsigned int y = cluster->pixels[j].y;
				const unsigned int z = cluster->pixels[j].z;
				const unsigned int index = (x + (y + dim[1] * z) * dim[0]) * 3;
				outPtr[index] = r;
				outPtr[index + 1] = g;
				outPtr[index + 2] = b;
			}
		}
	}
}

void vtkSuperpixelFilter::calcRandRgbExcludeZero(vtkImageData* output)
{
	float* outPtr = static_cast<float*>(output->GetScalarPointer());
	int* dim = output->GetDimensions();
	// For every cluster
	for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++)
	{
		Cluster* cluster = &clusters[i];
		// If it is still in the heap
		if (cluster->energy != -1.0f)
		{
			if (cluster->sumG > 0.0f)
			{
				// Create a random color and assign it to every pixel in the cluster
				float r = static_cast<float>(rand() % 255);
				float g = static_cast<float>(rand() % 255);
				float b = static_cast<float>(rand() % 255);
				for (unsigned int j = 0; j < cluster->pixels.size(); j++)
				{
					const unsigned int x = cluster->pixels[j].x;
					const unsigned int y = cluster->pixels[j].y;
					const unsigned int z = cluster->pixels[j].z;
					const unsigned int index = (x + (y + dim[1] * z) * dim[0]) * 3;
					outPtr[index] = r;
					outPtr[index + 1] = g;
					outPtr[index + 2] = b;
				}
			}
			else
			{
				// If it was never considered in the segmentation then just assign 0
				const unsigned int x = clusters->pixels[0].x;
				const unsigned int y = cluster->pixels[0].y;
				const unsigned int z = cluster->pixels[0].z;
				const unsigned int index = (x + (y + dim[1] * z) * dim[0]) * 3;
				outPtr[index] = 0.0f;
				outPtr[index + 1] = 0.0f;
				outPtr[index + 2] = 0.0f;
			}
		}
	}
}


unsigned int vtkSuperpixelFilter::initClusters(vtkImageData* input)
{
	float* inPtr = static_cast<float*>(input->GetScalarPointer());
	int* dim = input->GetDimensions();
	clusters = new Cluster[dim[0] * dim[1] * dim[2]];
	int index = 0;
	unsigned int zeroCount = 0;
	for (int z = 0; z < dim[2]; z++)
	{
		for (int y = 0; y < dim[1]; y++)
		{
			for (int x = 0; x < dim[0]; x++)
			{
				float g = inPtr[index];
				if (g <= 0.0f)
					zeroCount++;
				// Create a new node
				PixelNode node(x, y, z, g * colorWeight);
				// Each cluster is a collection of PixelNodes we start with one cluster for every PixelNode then merge clusters
				clusters[index] = Cluster();
				// Add the pixel to the cluster
				clusters[index].pixels.push_back(node);
				clusters[index].calcEnergy();
				index++;
			}
		}
	}

	return zeroCount;
}

// Adds the pairs to the min heap (if 1 is passed in for depth it will only setup 2d connections)
void vtkSuperpixelFilter::initHeap(vtkImageData* input)
{
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
					minHeap.insert(pair, -pair->dEnergy);
				}
			}
		}
	}
}

void vtkSuperpixelFilter::initHeapExcludeZero(vtkImageData* input)
{
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
				if (clusters[index].pixels[0].g > 0.0f && clusters[index1].pixels[0].g > 0.0f)
				{
					ClusterPair* pair = new ClusterPair(&clusters[index], &clusters[index1]);
					clusters[index].addEdge(pair);
					clusters[index1].addEdge(pair);
					pair->calcMergingCost();
					minHeap.insert(pair, -pair->dEnergy);
				}
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
				if (clusters[index].pixels[0].g > 0.0f && clusters[index1].pixels[0].g > 0.0f)
				{
					ClusterPair* pair = new ClusterPair(&clusters[index], &clusters[index1]);
					clusters[index].addEdge(pair);
					clusters[index1].addEdge(pair);
					pair->calcMergingCost();
					minHeap.insert(pair, -pair->dEnergy);
				}
			}
		}
	}

	// If depth is 1 then we are working with a 2d image and shouldn't add these pairs
	if (depth > 1)
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
					if (clusters[index].pixels[0].g > 0.0f && clusters[index1].pixels[0].g > 0.0f)
					{
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
}