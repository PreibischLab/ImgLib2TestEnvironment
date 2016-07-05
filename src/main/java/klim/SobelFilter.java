package klim;

import java.io.File;
import java.util.Iterator;

import javax.swing.text.View;

import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.RandomAccessiblePair;
import net.imglib2.view.Views;
import stephan.DistanceTransform;
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

	// set the values for the Sobel kernel
	// the values are set only for one axis
	// to get other kernels one should rotate 
	// the initial kernel 
	// the normalization is provided but 
	// it is not crucial
	public static void getKernelValues(int dim, int kType, float[] kValues){
		if (dim == 2){
			if (kType == 0){ // pi/4
				kValues[0] = 0;
				kValues[1] = 1; 
				kValues[2] = 2;
				kValues[3] = -1;
				kValues[4] = 0;
				kValues[5] = 1;
				kValues[6] = -2;
				kValues[7] = 1;
				kValues[8] = 0;
			}
			else{ // pi/2
				kValues[0] = 1;
				kValues[1] = 2; 
				kValues[2] = 1;
				kValues[3] = 0;
				kValues[4] = 0;
				kValues[5] = 0;
				kValues[6] = -1;
				kValues[7] = -2;
				kValues[8] = -1;
			}
			// normalization
			for (int i = 0; i < kValues.length; ++i) {
				kValues[i] /= 4;
			}		
		} 
		if (dim == 3){
			float[] kernel = new float[]{
					1,2,1, 		2,4,2, 		1,2,1,
					0,0,0, 		0,0,0,		0,0,0,
					-1,-2,-1, -2,-4,-2, -1,-2,-1};
			for (int i = 0; i < kValues.length; ++i) {
				kValues[i] = kernel[i];
				kValues[i] /= 16;
			}
		}
		// otherwise we don't know what to do
		// exit program
		if (dim != 2 && dim != 3){
			System.out.println("dimensionality is wrong");
			System.exit(1);
		}
		// return new float[]{};
	}

	public static  < T extends RealType< T >> void applyMedianFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2){
		final int n = img.numDimensions();
		final int[] mD = new int[n];
		// here we consider symmetric kernel
		for (int d = 0; d < n; ++d)
			mD[d] = 3;
		MedianFilter.medianFilter(img, img2, mD);
	}

	public static < T extends RealType< T >> void  applyGaussianFilter(RandomAccessibleInterval <T> img, RandomAccessibleInterval <T> img2) throws IncompatibleTypeException{
		double[] sigma = new double[ img.numDimensions() ];
		for (int d = 0; d < img.numDimensions(); ++d)
			sigma[d] = 4; // size of the radius
		Gauss3.gauss(sigma, Views.extendMirrorSingle(img), img2);
	}

	public static Img <FloatType> setKernel(int n){
		int kType = 0; // defines the direction of the sobel filter
		// fill in the kernel with proper values
		float[] kernelValues = new float [(int) Math.pow(3, n)]; 
		
		//System.out.println(kernelValues.length);
		
		getKernelValues(n, kType, kernelValues); 
	
		long [] kernelDimensions = new long[n];

		for (int d = 0; d < n; d++) 
			kernelDimensions[d] = 3; // the value is always set to 3 = size of the stencil
		
		// TODO: think which directions you need for a 3D case 

		// convert kernel to image 
		return ArrayImgs.floats( kernelValues, kernelDimensions);
	}

	// naive copy function
	public static <T extends RealType<T>> void copy(IterableInterval<T> in, RandomAccessibleInterval<T> out){
		Cursor<T> cursor = in.cursor();
		RandomAccess<T> randomAccess = out.randomAccess();
		
		while (cursor.hasNext()){
			cursor.fwd();
			randomAccess.setPosition(cursor);
			randomAccess.get().set(cursor.get().copy());		
		}
	}
	
	public static <T extends RealType<T>> void applySobelFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2, RandomAccessibleInterval<T> kernel, T minValue, T maxValue){
		RandomAccessibleInterval<T> tmp = img; // TODO: what is the right way to initialize it ?!
	
		// apply Sobel filter # of dimensions times
		for (int d = 0; d < img.numDimensions(); d++) {
			copy(Views.iterable(img), tmp);
			// convlove with each kernel
			new FFTConvolution<T>(tmp, Views.rotate(kernel, 0, d), new ArrayImgFactory<ComplexFloatType>()).convolve();

			// here we copy data to the destination image
			Cursor<T> tmpCursor = Views.iterable(tmp).cursor();
			RandomAccess<T> dstRandomAccess = img2.randomAccess();

			while (tmpCursor.hasNext()){
				tmpCursor.fwd();
				dstRandomAccess.setPosition(tmpCursor);
				T val = tmpCursor.get().copy();
				//System.out.println(val);
				val.mul(val);
				//System.out.println(val);
				dstRandomAccess.get().add(val);
			}

		}

		Cursor<T> dstCursor = Views.iterable(img2).cursor();
		while (dstCursor.hasNext()){
			dstCursor.fwd();	
			dstCursor.get().setReal(Math.sqrt(dstCursor.get().getRealFloat()));
		}
	
		Normalize.normalize(Views.iterable(img2), minValue, maxValue);

		// ImageJFunctions.show(img2);

	}

	public static <T extends RealType<T>> void distanceTransformKDTree(PointSampleList<T> worm, final RandomAccessibleInterval< BitType > img5, final RandomAccessibleInterval< FloatType > img6){

		final PointSampleList< BitType > oneList = new PointSampleList< BitType >( img5.numDimensions() );
		Cursor<T> wormCursor = worm.cursor(); 
		final RandomAccess<BitType> rImg5 = img5.randomAccess();

		while(wormCursor.hasNext()){
			wormCursor.fwd();
			rImg5.setPosition(wormCursor);
			if (rImg5.get().get()) {
				oneList.add(new Point(wormCursor), new BitType(true));
			}
		}

		final KDTree<BitType> tree = new KDTree<BitType>(oneList);
		final NearestNeighborSearchOnKDTree< BitType > search = new NearestNeighborSearchOnKDTree< BitType >(tree);


		final RandomAccess<FloatType> rImg6 = img6.randomAccess();
		wormCursor.reset();

		while(wormCursor.hasNext()){
			wormCursor.fwd();
			rImg5.setPosition(wormCursor);
			rImg6.setPosition(wormCursor);
			// System.out.println(wormCursor.get().getRealFloat());

			if(rImg5.get().getInteger() == 0){  
				search.search(rImg5);

				rImg6.get().setReal(search.getDistance());

				//	System.out.println("Hello!");
			}
			else{
				rImg6.get().setZero();
			}

			if(rImg5.get().getInteger() == 1){
				rImg6.get().set(80);
			}

		}

	}

	// process the worm and return the distance transform of it 
	// fix the declaration
	public static <T extends RealType<T>, U extends RealType<U>> void processWorm(RandomAccessibleInterval<T> initialImg, RandomAccessibleInterval<T> filterImg, RandomAccessibleInterval<T> edgeImg, RandomAccessibleInterval<BitType> thresholdImg, Img <T> distanceImg,
			T minValue, T maxValue,
			T tVal,
			U minVal, U maxVal
			) throws IncompatibleTypeException{		
		
		Normalize.normalize(Views.iterable(initialImg), minValue, maxValue);
		final int n = initialImg.numDimensions();
		// ImageJFunctions.show(initialImg);	

		// applyMedianFilter(initialImg, filterImg);
		applyGaussianFilter(initialImg, filterImg);	
		Img<T> kernel = (Img<T>) setKernel(n);
		applySobelFilter(filterImg, edgeImg, kernel, minValue, maxValue);
		ImageJFunctions.show(filterImg).setTitle("filterImg: in function");

		// check this thing below 
		// TODO: looks like this assignment is not working 
		// .copy() doesn't help! 
		thresholdImg = Thresholder.threshold((Img<T>)edgeImg, tVal, true, 1);
		ImageJFunctions.show(thresholdImg).setTitle("thresholdImg: in function");
//
//		// @TODO: make a function for this part 
//		// -------------------------------------
//		// -------------------------------------
//
//		Img< FloatType > img3 = imgFactory.create( img, new FloatType() ); 
//		Img< FloatType > img4 = imgFactory.create( img, new FloatType() );
//		Img< BitType > img5 = bitFactory.create( img, new BitType() );
//		Img< FloatType > img6 = imgFactory.create( img, new FloatType() );
//		Img< BitType > img7 = bitFactory.create( img, new BitType() );
//
//		final Cursor< BitType > first = dst.localizingCursor();
//		final RandomAccess< FloatType > second = img3.randomAccess();
//
//		while(first.hasNext()){
//			first.fwd();
//			second.setPosition(first);
//			second.get().set(first.get().getRealFloat()*255);
//		}
//
//		applySobelFilter(img3, img4, kernel);
//		ImageJFunctions.show(img4).setTitle("img4");
//		// -------------------------------------
//		// -------------------------------------
//		// @END_TODO ---------------------------
//		
//		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(dst, new IntType())); 
//		BoundingBox.setLabeling(dst, labeling);
//
//		// ImageJFunctions.show(labeling.getIndexImg());	
//		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
//
//		PointSampleList<T> worm = new PointSampleList<T>(img.numDimensions());
//		BoundingBox.setPointSampleList(labeling, (RandomAccessible<T>)img, worm);
//
//
//		img5 = Thresholder.threshold(img4, tVal, true, 1);
//
//		// @TODO: make a function for this part 
//		// ------------------------------------ 
//		// ------------------------------------
//
//
//		distanceTransformKDTree(worm, img5, img6);;
//
//		ImageJFunctions.show(img6).setTitle("I hope this looks fine");	
	}
	
	
	public static <T extends RealType<T>>void main(String[] args) throws IncompatibleTypeException{
		new ImageJ();
		File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");

		Img<FloatType> initialImg = ImgLib2Util.openAs32Bit(file);	
		ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
		ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
		
		Img< FloatType > filterImg = imgFactory.create( initialImg, new FloatType() );
		Img< BitType > thresholdImg = bitFactory.create( initialImg, new BitType() );
		Img< FloatType > edgeImg = imgFactory.create( initialImg, new FloatType() );
		Img< FloatType > distanceImg = imgFactory.create( initialImg, new FloatType() );
		
		processWorm(initialImg, filterImg, edgeImg, thresholdImg, distanceImg,
				new FloatType((float) 0), new FloatType((float) 255), new FloatType((float) 64),
				new BitType(false), new BitType(true));
		
		ImageJFunctions.show(initialImg).setTitle("initialImg");
		ImageJFunctions.show(filterImg).setTitle("filterImg");
		ImageJFunctions.show(edgeImg).setTitle("edgeImg");
		ImageJFunctions.show(thresholdImg).setTitle("thresholdImg");
		
		
		// TODO: Uncomment the one below to get everything back
		
		
		// temporary pic for calculations 
//		ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
//		Img< FloatType > img2 = imgFactory.create( img, new FloatType() );
//
//		// threshold output image 
//		ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
//		Img< BitType > dst = bitFactory.create( img, new BitType() );
//
//		FloatType minValue = new FloatType();
//		FloatType maxValue = new FloatType();
//		minValue.set(0);
//		maxValue.set(255);		
//		Normalize.normalize(img, minValue, maxValue);
//		final int n = img.numDimensions();
//		ImageJFunctions.show(img);	
//
//		img = applyGaussianFilter(img);	
//		applyMedianFilter(img, img2);
//		Img<FloatType> kernel = setKernel(n);
//		applySobelFilter(img, img2, kernel);
//
//		// check this thing below 
//		
//		FloatType tVal = new FloatType();
//		tVal.set((float) 64.0); // TODO: automate search of this value 
//		BitType minV = new BitType();
//		minV.setZero();
//		BitType maxV = new BitType();
//		maxV.setOne();
//		dst = Thresholder.threshold(img2, tVal, true, 1);
//		// ImageJFunctions.show(dst); 
//
//		// @TODO: make a function for this part 
//		// -------------------------------------
//		// -------------------------------------
//
//		Img< FloatType > img3 = imgFactory.create( img, new FloatType() ); 
//		Img< FloatType > img4 = imgFactory.create( img, new FloatType() );
//		Img< BitType > img5 = bitFactory.create( img, new BitType() );
//		Img< FloatType > img6 = imgFactory.create( img, new FloatType() );
//		Img< BitType > img7 = bitFactory.create( img, new BitType() );
//
//		final Cursor< BitType > first = dst.localizingCursor();
//		final RandomAccess< FloatType > second = img3.randomAccess();
//
//		while(first.hasNext()){
//			first.fwd();
//			second.setPosition(first);
//			second.get().set(first.get().getRealFloat()*255);
//		}
//
//		applySobelFilter(img3, img4, kernel);
//		ImageJFunctions.show(img4).setTitle("img4");
//		// -------------------------------------
//		// -------------------------------------
//		// @END_TODO ---------------------------
//		
//		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(dst, new IntType())); 
//		BoundingBox.setLabeling(dst, labeling);
//
//		// ImageJFunctions.show(labeling.getIndexImg());	
//		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
//
//		PointSampleList<T> worm = new PointSampleList<T>(img.numDimensions());
//		BoundingBox.setPointSampleList(labeling, (RandomAccessible<T>)img, worm);
//
//
//		img5 = Thresholder.threshold(img4, tVal, true, 1);
//
//		// @TODO: make a function for this part 
//		// ------------------------------------ 
//		// ------------------------------------
//
//
//		distanceTransformKDTree(worm, img5, img6);;
//
//		ImageJFunctions.show(img6).setTitle("I hope this looks fine");		 
//		// ------------------------------------
//		// ------------------------------------
//		// @END_TODO ---------------------------
//
//		file = new File("../Documents/Useful/initial_worms_pics/1003-red-one-2.tif");
//		Img<FloatType> img8 = ImgLib2Util.openAs32Bit(file);		
//		// temporary pic for calculations 
//		ImageJFunctions.show(img8).setTitle("straignt worm");	
//
//		/*
//		 * Here comes the same part as for the deformed worm
//		 * you should change the way the straight worm looks
//		 * */
//
//		Img< FloatType > img9 = imgFactory.create( img8, new FloatType() );
//
//		// threshold output image 
//		// ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
//		Img< BitType > dst1 = bitFactory.create( img8, new BitType() );
//
//		minValue = new FloatType();
//		maxValue = new FloatType();
//		minValue.set(0);
//		maxValue.set(255);		
//		Normalize.normalize(img8, minValue, maxValue);
//		// final int n = img.numDimensions();
//		ImageJFunctions.show(img8);	
//
//		img8 = applyGaussianFilter(img8);	
//		applyMedianFilter(img8, img9);
//		// Img<FloatType> kernel = setKernel(n);
//		applySobelFilter(img8, img9, kernel);
//		// ImageJFunctions.show(img8);
//		// ImageJFunctions.show(img9);
//		// check this thing below 
//
//		tVal = new FloatType();
//		tVal.set((float) 120.0);
//		minV = new BitType();
//		minV.setZero();
//		maxV = new BitType();
//		maxV.setOne();
//		dst1 = Thresholder.threshold(img9, tVal, true, 1);
//		ImageJFunctions.show(dst1);
//
//		// plain with distance transform here 
//
//		// @TODO: make a function for this part 
//		// -------------------------------------
//		// -------------------------------------
//
//		Img< FloatType > img10 = imgFactory.create( img8, new FloatType() ); 
//		Img< FloatType > img11 = imgFactory.create( img8, new FloatType() );
//		Img< BitType > img12 = bitFactory.create( img8, new BitType() );
//		Img< FloatType > img13 = imgFactory.create( img8, new FloatType() );
//		Img< BitType > img14 = bitFactory.create( img8, new BitType() );
//
//		final Cursor< BitType > first1 = dst1.localizingCursor();
//		final RandomAccess< FloatType > second1 = img10.randomAccess();
//
//		while(first1.hasNext()){
//			first1.fwd();
//			second1.setPosition(first1);
//			second1.get().set(first1.get().getRealFloat()*255);
//		}
//
//		applySobelFilter(img10, img11, kernel);
//		ImageJFunctions.show(img11).setTitle("img11");
//		// -------------------------------------
//		// -------------------------------------
//		// @END_TODO ---------------------------
//
//		final ImgLabeling<Integer, IntType> labeling1 = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(dst1, new IntType())); 
//		BoundingBox.setLabeling(dst1, labeling1);
//		// ImageJFunctions.show(labeling.getIndexImg());	
//		Cursor<IntType> cursor1 = Views.iterable(labeling1.getIndexImg()).cursor();
//
//		PointSampleList<T> worm1 = new PointSampleList<T>(img8.numDimensions());
//		BoundingBox.setPointSampleList(labeling1, (RandomAccessible<T>)img8, worm1);
//
//
//		img12 = Thresholder.threshold(img11, tVal, true, 1);
//
//		// @TODO: make a function for this part 
//		// ------------------------------------ 
//		// ------------------------------------
//
//
//		distanceTransformKDTree(worm1, img12, img13);
//
//
//		// img7 = Thresholder.threshold(img6, new FloatType((float)70.0), true, 1);
//
//		ImageJFunctions.show(img13).setTitle("I hope this looks fine");

		// ---------------------
		// ---------------------
		// ---------------------

		System.out.println("DONE!");
	}
}
