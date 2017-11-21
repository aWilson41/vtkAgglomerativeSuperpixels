#pragma once

#include "ClusterPair.h"
#include "Cluster.h"
#include "Mx\MxHeap.h"

#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vector>

class vtkSuperpixelFilter : public vtkImageAlgorithm
{
public:
	static vtkSuperpixelFilter* New();
	vtkTypeMacro(vtkSuperpixelFilter, vtkImageAlgorithm);

	vtkSuperpixelFilter() { }

	void SetInputData(vtkImageData* data) { vtkImageAlgorithm::SetInputData(data); }
	void SetNumberOfSuperpixels(unsigned int numSuperpixels) { vtkSuperpixelFilter::numSuperpixels = numSuperpixels; }

protected:
	int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

private: // Template code?
	vtkSuperpixelFilter(const vtkSuperpixelFilter&); // Not implemented
	void operator=(const vtkSuperpixelFilter&); // Not implemented

	void initClusters(unsigned char* inPtr, int width, int height, int depth);
	void initHeap(int width, int height, int depth, int numDim);
	void removeEdges(ClusterPair* pair, Cluster* c1, Cluster* c2);

private:
	// Min heap of cluster pairs each containing two clusters from the clusters array
	MxHeap minHeap;
	// Clusters array containg every cluster
	Cluster* clusters;
	unsigned int numSuperpixels = -1;
};