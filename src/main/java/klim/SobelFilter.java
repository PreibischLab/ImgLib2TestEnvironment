package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class SobelFilter {
	// this class applies Sobel Filter
	public static <T extends RealType<T>> void sobelFilter(
			final RandomAccessibleInterval< T > src, final RandomAccessibleInterval< T > dst, final int[] kernelDim){
		sobelFilter(Views.extendMirrorSingle(src), src, dst, kernelDim);
	}
	
	public static < T extends RealType<T>> void sobelFilter(
			final RandomAccessible< T > infSrc, final Interval srcInterval, final RandomAccessibleInterval< T > dst, final int[] kernelDim){
		
		final RandomAccessibleInterval<T> src = Views.interval(infSrc, srcInterval);
		
	}
	
	public static <T extends RealType<  T > > void getKernel(final int dim, final RandomAccessibleInterval<T> kernel){
		if (dim == 2){
			
		}
	}
	
	// set the values for the Sobel kernel
	// the values are set only for one axis
	// to get other kernels one should rotate 
	// the initial kernel 
	// the normalization is provided but 
	// it is not crucial
	public static float[] getKernelValues(int dim, int kType){
		float[] k = new float[9];
		if (dim == 2){
			if (kType == 0){
			k[0] = 0;
			k[1] = 1; 
			k[2] = 2;
			k[3] = -1;
			k[4] = 0;
			k[5] = 1;
			k[6] = -2;
			k[7] = 1;
			k[8] = 0;
		}
		else{
			k[0] = 1;
			k[1] = 2; 
			k[2] = 1;
			k[3] = 0;
			k[4] = 0;
			k[5] = 0;
			k[6] = -1;
			k[7] = -2;
			k[8] = -1;
		}

		for (int i = 0; i < k.length; ++i) {
				k[i] /= 4;
			}
			return k;			
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
		// new ImageJ();
		File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
		 Img<FloatType> img = ImgLib2Util.openAs32Bit(file);		
		// temporary pic for calculations 
		ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
		Img< FloatType > img2 = imgFactory.create( img, new FloatType() );
		
		// threshold output image 
		ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
		Img< BitType > dst = bitFactory.create( img, new BitType() );
		
		FloatType minValue = new FloatType();
		FloatType maxValue = new FloatType();
		minValue.set(0);
		maxValue.set(255);		
		Normalize.normalize(img, minValue, maxValue);
		
		double[] sigma = new double[ img.numDimensions() ];
		for (int d = 0; d < img.numDimensions(); ++d)
			sigma[d] = 4; // size of the radius
		img = Gauss.toFloat(sigma, img);	
		
		final int n = img.numDimensions();
		final int[] mD = new int[n];
		// here we consider symmetric kernel
		for (int d = 0; d < n; ++d)
			mD[d] = 3;
		MedianFilter.medianFilter(img, img2, mD);
		
		
		for (int kType = 0; kType <= 0; ++kType){

			// fill in the kernel with proper values
			float[] kernelValues = getKernelValues(img.numDimensions(), kType); 
			long [] kernelDimensions = new long[img.numDimensions()];

			for (int d = 0; d < img.numDimensions(); d++) 
				kernelDimensions[d] = 3; // the value is always set to 3 = size of the stencil


			// @TODO: think which directions you need for a 3D case 

			// convert kernel to image 
			Img<FloatType> kernel = ArrayImgs.floats( kernelValues, kernelDimensions);
			Img<FloatType> tmp = img.factory().create(img, new FloatType());	


			// apply Sobel filter # of dimensions times
			for (int d = 0; d < img.numDimensions(); d++) {
				tmp = img.copy();
				// convlove with each kernel
				new FFTConvolution<FloatType>(tmp, Views.rotate(kernel, 0, d), new ArrayImgFactory<ComplexFloatType>()).convolve();

				// Normalize.normalize(tmp, minValue, maxValue);
				// ImageJFunctions.show(tmp);

				// ImagePlus imgTmp = ImageJFunctions.show(tmp);
				// imgTmp.setDisplayRange(0, 255);
				// imgTmp.updateAndDraw();
				// enough to see the axis gradient

				// here we copy data to the destination image
				Cursor<FloatType> tmpCursor = tmp.cursor();
				RandomAccess<FloatType> dstRandomAccess = img2.randomAccess();

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
		}
		
		Cursor<FloatType> dstCursor = img2.cursor();
		while (dstCursor.hasNext()){
			dstCursor.fwd();	
			dstCursor.get().set((float) Math.sqrt(dstCursor.get().get() < 1e-6 ? 0 : dstCursor.get().get()));
			// dstCursor.get().set((float) Math.sqrt(dstCursor.get().get()));
			
		}
		
		Normalize.normalize(img2, minValue, maxValue);
		FloatType tVal = new FloatType();
		// tVal.set((float) 120);
		for (float thresholdVal = 55; thresholdVal <= 60; thresholdVal += 1){
			tVal.set((float) thresholdVal);
			BitType minV = new BitType();
			minV.setZero();
			BitType maxV = new BitType();
			maxV.setOne();
			dst = Thresholder.threshold(img2, tVal, true, 1);
			Normalize.normalize(dst, minV, maxV);
			ImageJFunctions.show(dst);
		}
		System.out.println("DONE!");
	}
}
