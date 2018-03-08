#pragma once
#include "Mx/MxHeap.h"

#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vector>

class Cluster;
class ClusterPair;

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

	void initClusters(vtkImageData* input);
	void initHeap(vtkImageData* input);
	void removeEdges(ClusterPair* pair);

	void calcColorLabels(vtkImageData* output);
	void calcRandRgb(vtkImageData* output);
	void calcAvgColors(vtkImageData* output);

private:
	// Min heap of cluster pairs each containing two clusters from the clusters array
	MxHeap minHeap;
	// Clusters array containing every cluster
	Cluster* clusters;
	unsigned int numSuperpixels = 0;
	float colorWeight = 1.0f;
	OutputType outputType = LABEL;
	bool swap = false;
};