package klim;

import java.io.File;

import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.filter.Rotator;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import stephan.MeanFilter;
import util.ImgLib2Util;

public class Sandbox {

	public Sandbox(){ 
		// set the file name 
		File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow.tif");
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");
		//File file = new File("src/main/resources/Bikesgray.jpg");
		
		// get + create initial ad final images 
		Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		Img <FloatType> dst = img.factory().create(img, img.firstElement());
		 
		// this part is necessary for normalization // 
		// TODO: Generic way? 
		FloatType minValue = new FloatType();
		FloatType maxValue = new FloatType();
		minValue.set(0);
		maxValue.set(255);		
		Normalize.normalize(img, minValue, maxValue);
		ImageJFunctions.show(img);
		// ---------------------------------------- // 
				
		// fill in the kernel with proper values
		float[] kernelValues = getKernelValues(img.numDimensions()); 
		long [] kernelDimensions = new long[img.numDimensions()];
		
		for (int d = 0; d < img.numDimensions(); d++) {
			kernelDimensions[d] = 3; // the value is always set to 3 = size of the stencil
		}
		// convert kernel to image 
		Img<FloatType> kernel = ArrayImgs.floats( kernelValues, kernelDimensions);
		Img<FloatType> tmp = img.factory().create(img, new FloatType());	
		// apply Sobel filter # of dimensions times
		for (int d = 0; d < img.numDimensions(); d++) {
			tmp = img.copy();
			// TODO: Check that the rotation is done correctly
			// TODO: Rotation is definitely wrong 
			// TODO: Maybe not the rotation but noise
			// TODO: Everything is fine 
			// the intensities map was wrong 
			// ---- 
			new FFTConvolution<FloatType>(tmp, Views.rotate(kernel, 0, d), new ArrayImgFactory<ComplexFloatType>()).convolve();
			
			
			// Normalize.normalize(tmp, minValue, maxValue);
			ImageJFunctions.show(tmp);
			
			
			// ImagePlus imgTmp = ImageJFunctions.show(tmp);
			// imgTmp.setDisplayRange(0, 255);
			// imgTmp.updateAndDraw();
			// enough to see the axis gradient
	
			// here we copy data to the destination image
			Cursor<FloatType> tmpCursor = tmp.cursor();
			RandomAccess<FloatType> dstRandomAccess = dst.randomAccess();
			
			while (tmpCursor.hasNext()){
				tmpCursor.fwd();
				dstRandomAccess.setPosition(tmpCursor);
				FloatType val = tmpCursor.get();
				//System.out.println(val);
				val.mul(val);
				//System.out.println(val);
				dstRandomAccess.get().add(val);
			}
					
		}
		
		
		// calculate the square root
		Cursor<FloatType> dstCursor = dst.cursor();
		while (dstCursor.hasNext()){
			dstCursor.fwd();
			
			
			// remove the noise from the picture 
			// naive way for that
			// dstCursor.get().set((float) Math.sqrt(dstCursor.get().get() < 1e-6 ? 0 : dstCursor.get().get()));
			dstCursor.get().set((float) Math.sqrt(dstCursor.get().get()));
			
		}
		ImageJFunctions.show(dst);		
	}
	
	
	// set the values for the Sobel kernel
	// the values are set only for one axis
	// to get other kernels one should rotate 
	// the initial kernel 
	// the normalization is provided but 
	// it is not crucial
	float[] getKernelValues(int dim){
		if (dim == 2){
			float[] kernel = new float[]{
					1,2,1, 0,0,0, -1,-2,-1};
			for (int i = 0; i < kernel.length; ++i) {
				kernel[i] /= 4;
			}
			return kernel;			
		} 
		if (dim == 3){
			float[] kernel = new float[]{
					1,2,1, 		2,4,2, 		1,2,1,
					0,0,0, 		0,0,0,		0,0,0,
					-1,-2,-1, -2,-4,-2, -1,-2,-1};
			for (int i = 0; i < kernel.length; ++i) {
				kernel[i] /= 16;
			}
		
			return kernel;
			
				
		}
		// otherwise we don't know what to do
		// exit program
		System.out.println("dimensionality is wrong");
		System.exit(1);
		return new float[]{};
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Sandbox();
		System.out.println("DONE!");
	}
}
