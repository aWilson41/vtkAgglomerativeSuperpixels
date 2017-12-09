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
	// If set to true zero values won't be included in merging and recieve a label of 0
	void SetExcludeZero(bool value) { excludeZero = value; }
	// If set to true it will output averaged color values instead of labels
	void SetOutputAvg(bool value) { outputAvg = value; }
	void SetNumberOfSuperpixels(unsigned int numSuperpixels) { vtkSuperpixelFilter::numSuperpixels = numSuperpixels; }
	// The color is scaled by this when added to the heap
	void SetWeight(double weight) { colorWeight = weight; }

protected:
	int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

private: // Template vtk filter code. why?
	vtkSuperpixelFilter(const vtkSuperpixelFilter&); // Not implemented
	void operator=(const vtkSuperpixelFilter&); // Not implemented

	unsigned int initClusters(float* inPtr, int width, int height, int depth = 1);
	void initHeap(int width, int height, int depth = 1);
	void initHeapExcludeZero(int width, int height, int depth = 1);
	void removeEdges(ClusterPair* pair, Cluster* c1, Cluster* c2);
	void calcColorLabels(float* outPtr, int width, int height, int depth = 1);
	void calcAvgColors(float* outPtr, int width, int height, int depth = 1);

private:
	// Min heap of cluster pairs each containing two clusters from the clusters array
	MxHeap minHeap;
	// Clusters array containing every cluster
	Cluster* clusters;
	unsigned int numSuperpixels = 0;
	double colorWeight = 1.0;
	// If this flag is turned on then 0's are ignored.
	bool excludeZero = false;
	bool outputAvg = false;
};