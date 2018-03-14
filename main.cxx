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
#include <vtkPNGWriter.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageLuminance.h>

// Slider
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSphereSource.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>

//void test3DImage();
void test2DImage();

vtkSmartPointer<vtkImageViewer2> imageViewer;
vtkSmartPointer<vtkSuperpixelFilter> superpixelFilter;
int numIter = 33;
int numSuperpixels = 500;

// Callback for slider for 3d image tests
//class vtkSliderCallback : public vtkCommand
//{
//public:
//	static vtkSliderCallback* New() { return new vtkSliderCallback; }
//	
//	virtual void Execute(vtkObject* caller, unsigned long, void*)
//	{
//		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
//		imageViewer->SetSlice(static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue());
//	}
//	vtkSliderCallback() : imageViewer(0) { }
//	vtkImageViewer2* imageViewer;
//};

// Define interaction style
//class vtkKeyPressInteractorStyle : public vtkInteractorStyleImage
//{
//public:
//	static vtkKeyPressInteractorStyle* New();
//	vtkTypeMacro(vtkKeyPressInteractorStyle, vtkInteractorStyleImage);
//
//	virtual void OnKeyPress()
//	{
//		// Get the keypress
//		vtkRenderWindowInteractor* rwi = this->Interactor;
//		std::string key = rwi->GetKeySym();
//
//		// Handle an arrow key
//		if (key == "space")
//		{
//			superpixelFilter->SetSwapIterations(++numIter);
//			/*numSuperpixels /= 2;
//			superpixelFilter->SetNumberOfSuperpixels(numSuperpixels);*/
//			superpixelFilter->Update();
//
//			// To write the image as png we cast to uchar
//			vtkSmartPointer<vtkImageCast> writeCast = vtkSmartPointer<vtkImageCast>::New();
//			writeCast->SetInputData(superpixelFilter->GetOutput());
//			writeCast->SetOutputScalarTypeToUnsignedChar();
//			writeCast->Update();
//			vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
//			writer->SetInputData(writeCast->GetOutput());
//			writer->SetFileName(("C:/Users/Andx_/Desktop/spOutput/spOutput" + std::to_string(numIter) + ".png").c_str());
//			writer->Write();
//
//			imageViewer->Render();
//			printf("Update\n");
//		}
//
//		// Forward events
//		vtkInteractorStyleTrackballCamera::OnKeyPress();
//	}
//
//};
//vtkStandardNewMacro(vtkKeyPressInteractorStyle);

int main(int argc, char* argv[])
{
	test2DImage();
	//test3DImage();

	return EXIT_SUCCESS;
}

void test2DImage()
{
	// Read the 2d png
	vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
	reader->SetFileName("C:/Users/Andx_/Desktop/8049.png");
	reader->Update();
	// Cast to float in case it's not already float (since my filter only works with float)
	vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
	cast->SetInputData(reader->GetOutput());
	cast->SetOutputScalarTypeToFloat();
	cast->Update();

	vtkSmartPointer<vtkImageData> input = cast->GetOutput();
	if (input->GetNumberOfScalarComponents() > 1)
	{
		vtkSmartPointer<vtkImageLuminance> toGrayScale = vtkSmartPointer<vtkImageLuminance>::New();
		toGrayScale->SetInputData(cast->GetOutput());
		toGrayScale->Update();
		input = toGrayScale->GetOutput();
	}

	// Superpixel segment
	superpixelFilter = vtkSmartPointer<vtkSuperpixelFilter>::New();
	superpixelFilter->SetInputData(input);
	superpixelFilter->SetNumberOfSuperpixels(numSuperpixels);
	superpixelFilter->SetSwap(true);
	superpixelFilter->SetSwapIterations(numIter);
	superpixelFilter->SetOutputType(vtkSuperpixelFilter::AVGCOLOR);
	superpixelFilter->Update();

	// Visualize
	imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	imageViewer->GetRenderWindow()->SetSize(1000, 800);
	imageViewer->SetInputData(superpixelFilter->GetOutput());
	imageViewer->GetImageActor()->InterpolateOff();
	imageViewer->SetupInteractor(renderWindowInteractor);
	//renderWindowInteractor->SetInteractorStyle(vtkSmartPointer<vtkKeyPressInteractorStyle>::New());

	// Render image
	imageViewer->Render();
	renderWindowInteractor->Start();
}

//void test3DImage()
//{
//	// Read 2d or 3d meta image (mhd)
//	vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
//	reader->SetFileName("MPRAGE_iso_sag.mhd");
//	//reader->SetFileName("C:/Users/Andx_/Desktop/3DisoFLAIR_sagXY.mhd");
//	//reader->SetFileName("C:/Users/Andx_/Desktop/ROI/3DisoFLAIR_sagGaussian.mhd");
//	reader->Update();
//
//	// Grab the first component forcing rgb/lab images to grayscale
//	vtkSmartPointer<vtkImageExtractComponents> extractComp = vtkSmartPointer<vtkImageExtractComponents>::New();
//	extractComp->SetInputData(reader->GetOutput());
//	extractComp->SetComponents(0);
//	extractComp->Update();
//
//	// Cast to float in case it's not already float (since my filter only works with float)
//	vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
//	cast->SetInputData(extractComp->GetOutput());
//	cast->SetOutputScalarTypeToFloat();
//	cast->Update();
//
//	// Superpixel segment
//	vtkSmartPointer<vtkSuperpixelFilter> superpixelFilter = vtkSmartPointer<vtkSuperpixelFilter>::New();
//	superpixelFilter->SetInputData(cast->GetOutput());
//	superpixelFilter->SetNumberOfSuperpixels(1000);
//	superpixelFilter->SetOutputType(vtkSuperpixelFilter::RANDRGB);
//	superpixelFilter->Update();
//
//	vtkSmartPointer<vtkImageCast> imageCast = vtkSmartPointer<vtkImageCast>::New();
//	imageCast->SetInputData(superpixelFilter->GetOutput());
//	imageCast->SetOutputScalarTypeToUnsignedChar();
//	imageCast->Update();
//
//	// Visualize
//	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	imageViewer->GetRenderWindow()->SetSize(1000, 800);
//	imageViewer->SetInputData(imageCast->GetOutput());
//	imageViewer->SetSlice((imageViewer->GetSliceMax() + imageViewer->GetSliceMin()) / 2); // Set to middle slice
//	imageViewer->GetImageActor()->InterpolateOff();
//	imageViewer->SetupInteractor(renderWindowInteractor);
//
//#pragma region Setup Slider
//	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
//	sliderRep->SetMinimumValue(0);
//	sliderRep->SetMaximumValue(superpixelFilter->GetOutput()->GetDimensions()[2]);
//	sliderRep->SetValue(90);
//	sliderRep->GetSliderProperty()->SetColor(1, 0, 0);
//	sliderRep->GetLabelProperty()->SetColor(1, 0, 0);
//	sliderRep->GetSelectedProperty()->SetColor(0, 1, 0);
//	sliderRep->GetTubeProperty()->SetColor(1, 1, 0);
//	sliderRep->GetCapProperty()->SetColor(1, 1, 0);
//	sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
//	sliderRep->GetPoint1Coordinate()->SetValue(40, 40);
//	sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
//	sliderRep->GetPoint2Coordinate()->SetValue(200, 40);
//	vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
//	sliderWidget->SetInteractor(renderWindowInteractor);
//	sliderWidget->SetRepresentation(sliderRep);
//	sliderWidget->SetAnimationModeToAnimate();
//	sliderWidget->EnabledOn();
//
//	vtkSmartPointer<vtkSliderCallback> callback = vtkSmartPointer<vtkSliderCallback>::New();
//	callback->imageViewer = imageViewer;
//	sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);
//#pragma endregion
//
//	// Render image
//	imageViewer->Render();
//	/*imageViewer->GetRenderer()->ResetCamera();
//	imageViewer->Render();*/
//	renderWindowInteractor->Start();
//}