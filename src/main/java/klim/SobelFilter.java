package klim;

import java.io.File;
import java.util.Iterator;

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
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
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

	public static  < T extends RealType< T >> void applyMedianFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2){
		final int n = img.numDimensions();
		final int[] mD = new int[n];
		// here we consider symmetric kernel
		for (int d = 0; d < n; ++d)
			mD[d] = 3;
		MedianFilter.medianFilter(img, img2, mD);
	}

	public static Img <FloatType> applyGaussianFilter(Img <FloatType> img){
		double[] sigma = new double[ img.numDimensions() ];
		for (int d = 0; d < img.numDimensions(); ++d)
			sigma[d] = 4; // size of the radius
		// Gauss3.gauss(sigma, source, target);
		return Gauss.toFloat(sigma, img);
	}

	public static Img<FloatType> setKernel(int n ){
		int kType = 1; // defines the direction of the sobel filter
		// fill in the kernel with proper values
		float[] kernelValues = getKernelValues(n, kType); 
		long [] kernelDimensions = new long[n];

		for (int d = 0; d < n; d++) 
			kernelDimensions[d] = 3; // the value is always set to 3 = size of the stencil


		// @TODO: think which directions you need for a 3D case 

		// convert kernel to image 
		Img<FloatType> kernel = ArrayImgs.floats( kernelValues, kernelDimensions);
		return kernel;
	}

	public static void applySobelFilter(Img <FloatType> img, Img <FloatType> img2, Img<FloatType> kernel){
		Img<FloatType> tmp = img.factory().create(img, new FloatType());
		// apply Sobel filter # of dimensions times
		for (int d = 0; d < img.numDimensions(); d++) {
			tmp = img.copy();
			// convlove with each kernel
			new FFTConvolution<FloatType>(tmp, Views.rotate(kernel, 0, d), new ArrayImgFactory<ComplexFloatType>()).convolve();

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

		Cursor<FloatType> dstCursor = img2.cursor();
		while (dstCursor.hasNext()){
			dstCursor.fwd();	
			// dstCursor.get().set((float) Math.sqrt(dstCursor.get().get() < 1e-6 ? 0 : dstCursor.get().get()));
			dstCursor.get().set((float) Math.sqrt(dstCursor.get().get()));

		}
		FloatType minValue = new FloatType();
		FloatType maxValue = new FloatType();
		minValue.set(0);
		maxValue.set(255);		
		Normalize.normalize(img2, minValue, maxValue);

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

	public static <T extends RealType<T>>void main(String[] args){
		new ImageJ();
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
		final int n = img.numDimensions();
		ImageJFunctions.show(img);	

		img = applyGaussianFilter(img);	
		applyMedianFilter(img, img2);
		Img<FloatType> kernel = setKernel(n);
		applySobelFilter(img, img2, kernel);
		// ImageJFunctions.show(img2)

		// check this thing below 

		FloatType tVal = new FloatType();
		tVal.set((float) 64.0);
		BitType minV = new BitType();
		minV.setZero();
		BitType maxV = new BitType();
		maxV.setOne();
		dst = Thresholder.threshold(img2, tVal, true, 1);
		// ImageJFunctions.show(dst);

		// plain with distance transform here 

		// @TODO: make a function for this part 
		// -------------------------------------
		// -------------------------------------

		Img< FloatType > img3 = imgFactory.create( img, new FloatType() ); 
		Img< FloatType > img4 = imgFactory.create( img, new FloatType() );
		Img< BitType > img5 = bitFactory.create( img, new BitType() );
		Img< FloatType > img6 = imgFactory.create( img, new FloatType() );
		Img< BitType > img7 = bitFactory.create( img, new BitType() );

		final Cursor< BitType > first = dst.localizingCursor();
		final RandomAccess< FloatType > second = img3.randomAccess();

		while(first.hasNext()){
			first.fwd();
			second.setPosition(first);
			second.get().set(first.get().getRealFloat()*255);
		}

		applySobelFilter(img3, img4, kernel);
		ImageJFunctions.show(img4).setTitle("img4");
		// -------------------------------------
		// -------------------------------------
		// @END_TODO ---------------------------


		final ImgLabeling<Integer, IntType> labeling = BoundingBox.setLabeling(dst);
		// ImageJFunctions.show(labeling.getIndexImg());	
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();

		PointSampleList<T> worm = new PointSampleList<T>(img.numDimensions());
		BoundingBox.setPointSampleList(labeling, (RandomAccessible<T>)img, worm);


		img5 = Thresholder.threshold(img4, tVal, true, 1);

		// @TODO: make a function for this part 
		// ------------------------------------ 
		// ------------------------------------


		distanceTransformKDTree(worm, img5, img6);


		// img7 = Thresholder.threshold(img6, new FloatType((float)70.0), true, 1);

		ImageJFunctions.show(img6).setTitle("I hope this looks fine");	
		// ImageJFunctions.show(img7).setTitle("Here we go");	 
		// ------------------------------------
		// ------------------------------------
		// @END_TODO ---------------------------

		file = new File("../Documents/Useful/initial_worms_pics/1003-red-one-2.tif");
		Img<FloatType> img8 = ImgLib2Util.openAs32Bit(file);		
		// temporary pic for calculations 
		ImageJFunctions.show(img8).setTitle("straignt worm");	

		/*
		 * Here comes the same part as for the deformed worm
		 * you should change the way the straight worm looks
		 * */

		Img< FloatType > img9 = imgFactory.create( img8, new FloatType() );

		// threshold output image 
		// ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();
		Img< BitType > dst1 = bitFactory.create( img8, new BitType() );

		minValue = new FloatType();
		maxValue = new FloatType();
		minValue.set(0);
		maxValue.set(255);		
		Normalize.normalize(img8, minValue, maxValue);
		// final int n = img.numDimensions();
		ImageJFunctions.show(img8);	

		img8 = applyGaussianFilter(img8);	
		applyMedianFilter(img8, img9);
		// Img<FloatType> kernel = setKernel(n);
		applySobelFilter(img8, img9, kernel);
		// ImageJFunctions.show(img8);
		// ImageJFunctions.show(img9);
		// check this thing below 

		tVal = new FloatType();
		tVal.set((float) 120.0);
		minV = new BitType();
		minV.setZero();
		maxV = new BitType();
		maxV.setOne();
		dst1 = Thresholder.threshold(img9, tVal, true, 1);
		ImageJFunctions.show(dst1);

		// plain with distance transform here 

		// @TODO: make a function for this part 
		// -------------------------------------
		// -------------------------------------

		Img< FloatType > img10 = imgFactory.create( img8, new FloatType() ); 
		Img< FloatType > img11 = imgFactory.create( img8, new FloatType() );
		Img< BitType > img12 = bitFactory.create( img8, new BitType() );
		Img< FloatType > img13 = imgFactory.create( img8, new FloatType() );
		Img< BitType > img14 = bitFactory.create( img8, new BitType() );

		final Cursor< BitType > first1 = dst1.localizingCursor();
		final RandomAccess< FloatType > second1 = img10.randomAccess();

		while(first1.hasNext()){
			first1.fwd();
			second1.setPosition(first1);
			second1.get().set(first1.get().getRealFloat()*255);
		}

		applySobelFilter(img10, img11, kernel);
		ImageJFunctions.show(img11).setTitle("img11");
		// -------------------------------------
		// -------------------------------------
		// @END_TODO ---------------------------


		final ImgLabeling<Integer, IntType> labeling1 = BoundingBox.setLabeling(dst1);
		// ImageJFunctions.show(labeling.getIndexImg());	
		Cursor<IntType> cursor1 = Views.iterable(labeling1.getIndexImg()).cursor();

		PointSampleList<T> worm1 = new PointSampleList<T>(img8.numDimensions());
		BoundingBox.setPointSampleList(labeling1, (RandomAccessible<T>)img8, worm1);


		img12 = Thresholder.threshold(img11, tVal, true, 1);

		// @TODO: make a function for this part 
		// ------------------------------------ 
		// ------------------------------------


		distanceTransformKDTree(worm1, img12, img13);


		// img7 = Thresholder.threshold(img6, new FloatType((float)70.0), true, 1);

		ImageJFunctions.show(img13).setTitle("I hope this looks fine");

		// ---------------------
		// ---------------------
		// ---------------------

		// debug
		//		Img <FloatType> output = img.factory().create(img, img.firstElement());
		//		RandomAccess <FloatType> randomAccessOut = output.randomAccess();
		//		
		//		Cursor<T> it = worm.cursor(); 
		//		long [] pos = new long[2]; // change 2 to dimension of the image 
		//		while(it.hasNext()){
		//			it.fwd();
		//			// it.localize(pos);
		//			randomAccessOut.setPosition(it);
		//			randomAccessOut.get().set(it.get().getRealFloat());			
		//		}
		//		
		//		ImageJFunctions.show(output).setTitle("output");



		//		// ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
		//		Img< FloatType > imgBound = imgFactory.create( img, new FloatType() );
		//		Img< FloatType > imgBound2 = imgFactory.create( img, new FloatType() );
		//		// DistanceTransform.distanceTransformKD(dst, imgBound);
		//		
		//		final Cursor< BitType > first = dst.localizingCursor();
		//		final RandomAccess< FloatType > second = imgBound.randomAccess();
		//		
		//		
		//		while(first.hasNext()){
		//			first.fwd();
		//			second.setPosition(first);
		//			second.get().set(first.get().getRealFloat()*255);
		//		}
		//		
		//		// ImageJFunctions.show(imgBound);
		//		applySobelFilter(imgBound, imgBound2, kernel);
		//		// maxValue.set(10);
		//		// Normalize.normalize(imgBound, minValue, maxValue);
		//		ImageJFunctions.show(imgBound2);
		//		
		//		final Cursor< FloatType > third = imgBound2.localizingCursor();
		//		Img< BitType > imgBoundBit = bitFactory.create( img, new BitType() );
		//		final RandomAccess< BitType> forth = imgBoundBit.randomAccess();
		//		
		//		while(third.hasNext()){
		//			third.fwd();
		//			forth.setPosition(third);
		//			forth.get().set(third.get().getRealFloat() > 2.0 ? true : false);
		//		}
		//		
		//		ImageJFunctions.show(imgBoundBit);
		//		
		//		
		//		// DistanceTransform.distanceTransformKD(imgBoundBit, imgBound);
		//		ImageJFunctions.show(imgBound);




		System.out.println("DONE!");
	}
}
