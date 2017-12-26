#pragma once

#include "ClusterPair.h"
#include "Cluster.h"
#include "Mx/MxHeap.h"

#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vector>

class vtkSuperpixelFilter : public vtkImageAlgorithm
{
public:
	// Color labels: Sequential unique grayscale for every cluster (ie: 0, 1, 2, ...)
	// Random rgb: Random rgb value for every cluster
	// Average color: Averages the grayscale output
	enum OutputType
	{
		LABEL,
		RANDRGB,
		AVGCOLOR
	};

public:
	static vtkSuperpixelFilter* New();
	vtkTypeMacro(vtkSuperpixelFilter, vtkImageAlgorithm);

	vtkSuperpixelFilter() { }

	void SetInputData(vtkImageData* data) { vtkImageAlgorithm::SetInputData(data); }
	// If set to true zero values won't be included in merging and recieve a label of 0
	void SetExcludeZero(bool value) { excludeZero = value; }
	void SetOutputType(OutputType outputType) { vtkSuperpixelFilter::outputType = outputType; }
	void SetNumberOfSuperpixels(unsigned int numSuperpixels) { vtkSuperpixelFilter::numSuperpixels = numSuperpixels; }
	// The color is scaled by this when added to the heap
	void SetWeight(double weight) { colorWeight = weight; }
	// If set to true the output will be swapped
	void SetSwap(bool value) { swap = value; }

protected:
	int RequestInformation(vtkInformation* request, vtkInformationVector** inputVec, vtkInformationVector* outputVec) VTK_OVERRIDE;
	int RequestData(vtkInformation* request, vtkInformationVector** inputVec, vtkInformationVector* outputVec) VTK_OVERRIDE;

private: // Template vtk filter code. why?
	vtkSuperpixelFilter(const vtkSuperpixelFilter&); // Not implemented
	void operator=(const vtkSuperpixelFilter&); // Not implemented

	unsigned int initClusters(float* inPtr, int width, int height, int depth = 1);
	void initHeap(int width, int height, int depth = 1);
	void initHeapExcludeZero(int width, int height, int depth = 1);
	void removeEdges(ClusterPair* pair);
	void calcColorLabels(float* outPtr, int width, int height, int depth = 1);
	void calcRandRgb(float* outPtr, int width, int height, int depth = 1);
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
	OutputType outputType = LABEL;
	bool swap = false;
};