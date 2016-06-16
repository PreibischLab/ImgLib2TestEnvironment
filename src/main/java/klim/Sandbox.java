package klim;

import java.io.File;

import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.filter.Rotator;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.algorithm.fft.InverseFourierConvolution;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
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
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import stephan.MeanFilter;
import util.ImgLib2Util;

public class Sandbox {
	
	private static final boolean filterType = false;

	public Sandbox(){ 
		// set the file name 
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");
		File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
		// File file = new File("src/main/resources/Bikesgray.jpg");
		
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
		
		// preprocessing
		
		if (filterType){
			// gauss-smooth the picture a little bit
			double[] sigma = new double[ img.numDimensions() ];
			for (int d = 0; d < img.numDimensions(); ++d)
				sigma[d] = 10; // size of the radius
			 img = Gauss.toFloat(sigma, img);			
		}
		else{
			final int n = img.numDimensions();
			final long[] min = new long[n];
			final long[] max = new long[n];
	
			long medianSize = 1;
			for (int d = 0; d < n; d++) {
				min[d] = -medianSize;
				max[d] = medianSize;
			}
			int mD = 3;
			MedianFilter.medianFilter(img, dst, new int[]{mD, mD});
		}
		
		ImageJFunctions.show(img);
		
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
		
		ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
		Img< BitType > display = bitFactory.create( dst, new BitType() );
		
		FloatType tVal = new FloatType();
		if (filterType){
			tVal.set((float)0.8); // 
		
			
		}
		else{
			for(int i = 20; i < 40; i+= 2){
				tVal.set((float)(7.8 + i)); // 
				display = Thresholder.threshold(dst, tVal, true, 1);
				// ImageJFunctions.show(dst);	
				ImageJFunctions.show(display);	
			}

		}
		// display = Thresholder.threshold(dst, tVal, true, 1);
		// ImageJFunctions.show(dst);	
		// ImageJFunctions.show(display);		
		
		// Sharpen image
		// float[] kernelValuesS = new float[]{-1,-1,-1,
		//									-1, 25,-1,
		// 									-1,-1,-1};
		// Img<FloatType> kernelS = ArrayImgs.floats( kernelValuesS, kernelDimensions);
		// new FFTConvolution<FloatType>(dst, kernelS, new ArrayImgFactory<ComplexFloatType>()).convolve();
		// ImageJFunctions.show(dst);	
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
					0, 1, 2, -1, 0, 1, -2, 1, 0};
			
//			float[] kernel = new float[]{
//					1,2,1, 0,0,0, -1,-2,-1};
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
	
	public static <T extends Comparable<T>, U extends RealType<U>> Img <U> 
	findAndDisplaylocalMaxima(RandomAccessibleInterval<T> src, ImgFactory<U> imageFactory, U outputType){
		Img <U> output = imageFactory.create(src, outputType);
		Interval interval = Intervals.expand(src, -1);
		src = Views.interval(src, interval);
		final Cursor<T> center = Views.iterable(src).cursor();
		final RectangleShape shape = new RectangleShape(1, true);
		for (final Neighborhood<T> localNeighborhood: shape.neighborhoods(src)){
			final T centerValue = center.next();
			boolean isMax = true; 
			for (final T value : localNeighborhood){
				if (centerValue.compareTo(value) <= 0){
					isMax = false;
					break;
				}
			}
			if (isMax){
				HyperSphere<U> hyperSphere = new HyperSphere <U> (output, center, 1);
				for (U value : hyperSphere)
					value.setOne();
			}
		}
		return output;
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Sandbox();
		System.out.println("DONE!");
	}
}
