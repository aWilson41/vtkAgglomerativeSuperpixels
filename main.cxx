#include "vtkSuperpixelFilter.h"

#include <vtkSmartPointer.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageViewer2.h>
#include <vtkPNGReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageCast.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>

// Slider
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSphereSource.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>

class vtkSliderCallback : public vtkCommand
{
public:
	static vtkSliderCallback* New() { return new vtkSliderCallback; }
	
	virtual void Execute(vtkObject* caller, unsigned long, void*)
	{
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		imageViewer->SetSlice(static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue());
	}
	vtkSliderCallback() : imageViewer(0) { }
	vtkImageViewer2* imageViewer;
};

int main(int argc, char* argv[])
{
	// Read the 2d png
	vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
	//reader->SetFileName("C:/Users/Andx_/Desktop/testing.png");
	reader->SetFileName("test.png");
	reader->Update();

	// Read 2d or 3d meta image (mhd)
	//vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
	//reader->SetFileName("C:/Users/Andx_/Desktop/ROI/3DisoFLAIR_sagXYXYThresholded.mhd");
	////reader->SetFileName("C:/Users/Andx_/Desktop/3DisoFLAIR_sagXY.mhd");
	////reader->SetFileName("C:/Users/Andx_/Desktop/ROI/3DisoFLAIR_sagGaussian.mhd");
	//reader->Update();

	// Grab the first component forcing rgb/lab images to grayscale
	vtkSmartPointer<vtkImageExtractComponents> extractComp = vtkSmartPointer<vtkImageExtractComponents>::New();
	extractComp->SetInputData(reader->GetOutput());
	extractComp->SetComponents(0);
	extractComp->Update();

	// Cast to uchar in case it's not already uchar
	vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
	cast->SetInputData(extractComp->GetOutput());
	cast->SetOutputScalarTypeToUnsignedChar();
	cast->Update();

	// Superpixel segment
	vtkSmartPointer<vtkSuperpixelFilter> superpixelFilter = vtkSmartPointer<vtkSuperpixelFilter>::New();
	superpixelFilter->SetInputData(cast->GetOutput());
	superpixelFilter->SetNumberOfSuperpixels(100);
	superpixelFilter->Update();

	// Visualize
	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	imageViewer->GetRenderWindow()->SetSize(1000, 800);
	imageViewer->SetInputData(superpixelFilter->GetOutput());
	imageViewer->SetSlice(90);
	imageViewer->GetImageActor()->InterpolateOff();
	imageViewer->SetupInteractor(renderWindowInteractor);

	#pragma region Setup Slider
	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
	sliderRep->SetMinimumValue(0);
	sliderRep->SetMaximumValue(superpixelFilter->GetOutput()->GetDimensions()[2]);
	sliderRep->SetValue(90);
	sliderRep->GetSliderProperty()->SetColor(1, 0, 0);
	sliderRep->GetLabelProperty()->SetColor(1, 0, 0);
	sliderRep->GetSelectedProperty()->SetColor(0, 1, 0);
	sliderRep->GetTubeProperty()->SetColor(1, 1, 0);
	sliderRep->GetCapProperty()->SetColor(1, 1, 0);
	sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
	sliderRep->GetPoint1Coordinate()->SetValue(40, 40);
	sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
	sliderRep->GetPoint2Coordinate()->SetValue(200, 40);
	vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
	sliderWidget->SetInteractor(renderWindowInteractor);
	sliderWidget->SetRepresentation(sliderRep);
	sliderWidget->SetAnimationModeToAnimate();
	sliderWidget->EnabledOn();

	vtkSmartPointer<vtkSliderCallback> callback = vtkSmartPointer<vtkSliderCallback>::New();
	callback->imageViewer = imageViewer;
	sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);
	#pragma endregion

	// Render image
	imageViewer->Render();
	imageViewer->GetRenderer()->ResetCamera();
	imageViewer->Render();
	renderWindowInteractor->Start();

	// Write results
	vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
	writer->SetInputData(superpixelFilter->GetOutput());
	writer->SetFileName("output.mhd");
	writer->SetRAWFileName("outputRAW.raw");
	writer->Write();

	return EXIT_SUCCESS;
}