#pragma once
#include "Mx/MxHeap.h"

#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vector>

class Cluster;
class ClusterPair;
class PixelNode;

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
		AVGCOLOR,
		MAXCOLOR,
		MINCOLOR
	};

public:
	static vtkSuperpixelFilter* New();
	vtkTypeMacro(vtkSuperpixelFilter, vtkImageAlgorithm);

	vtkSuperpixelFilter() { }

	void SetInputData(vtkImageData* data) { vtkImageAlgorithm::SetInputData(data); }
	void SetOutputType(OutputType outputType) { vtkSuperpixelFilter::outputType = outputType; this->Modified(); }
	vtkSetMacro(NumberOfSuperpixels, unsigned int);
	vtkGetMacro(NumberOfSuperpixels, unsigned int);
	vtkSetMacro(SwapIterations, unsigned int);
	vtkGetMacro(SwapIterations, unsigned int);
	// The color is scaled by this when added to the heap
	vtkSetMacro(ColorWeight, double);
	vtkGetMacro(ColorWeight, double);

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
	void calcMaxColors(vtkImageData* output);
	void calcMinColors(vtkImageData* output);

	// Returns the largest change in energy of the swaps made
	void computeSwap(int width, int height, int depth);

private:
	// Min heap of cluster pairs each containing two clusters from the clusters array
	MxHeap* minHeap;
	// Clusters array containing every cluster made
	Cluster* clusters;
	// Points to all the resulting clusters (subset of clusters)
	std::vector<Cluster*> finalClusters;
	PixelNode* px;

	OutputType outputType = LABEL;
	unsigned int NumberOfSuperpixels = 0;
	double ColorWeight = 1.0f;
	unsigned int SwapIterations = 0;
};