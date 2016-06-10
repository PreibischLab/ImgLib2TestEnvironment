package klim;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import ij.ImageJ;

import ij.io.ImageWriter;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.MixedTransformView;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MedianFilter {
	
	
	private static final boolean debug = true;
	
	
	// TODO: write description 
	/*
	 * 
	 * */
	// TODO:  /*correct type declaration */
	public static < T extends RealType<  T > & Comparable<T> > void medianFilter(
			final RandomAccessibleInterval< T > src, final RandomAccessibleInterval< T > dst, final Interval kernel){

		final RandomAccessible<T> infSrc = Views.extendMirrorSingle(src);
		final Cursor<T> cSrc = Views.iterable(src).localizingCursor();
		final RandomAccess<T> rDst = dst.randomAccess();
		
		final int n = src.numDimensions();
		// store kernel boundaries
		final long[] min = new long[n];
		final long[] max = new long[n];
		for (int d = 0; d < kernel.numDimensions(); ++d){
			min[d] = cSrc.getLongPosition(d) - kernel.dimension(d)/2;
			max[d] = cSrc.getLongPosition(d) + kernel.dimension(d)/2;
		}
		// store previous/current position of cursor
		long[] pPos = new long[n]; 
		long[] cPos = new long[n]; 
		
		// contains all elements of the kernel
		IterableInterval<T> histogram = Views.interval(infSrc, min, max);
		// List<T> histogramList = new ArrayList<T>(); 
		List<T> histogramList = new LinkedList<T>(); 
		
		for (T h : histogram)
			histogramList.add(h.copy()); 					
		Collections.sort(histogramList);	
		
		// @DEBUG: delete idx after done 
		long idx = 1;
		while(cSrc.hasNext()){
			cSrc.localize(pPos);
			cSrc.fwd();
			cSrc.localize(cPos);
				
			// @TODO: Dirty! Clean this
			long[] tmp = new long[2];
			checkDist(pPos, cPos, tmp);	 // check if the cursor moved only by one step
			int direction = (int)tmp[0]; // -2 - too far; -1 - init value; >=0 - dimension for movement
			long step = tmp[1]; // shows the direction of movement (+1/-1)
			
			// define new boundaries of the new kernel-window
			for (int d = 0; d < n; ++d){
				min[d] = cSrc.getLongPosition(d) - kernel.dimension(d)/2;
				max[d] = cSrc.getLongPosition(d) + kernel.dimension(d)/2;
			}
			
			// adjust the histogram window
			histogram = Views.interval(infSrc, min, max);
						
			if (direction == -2){ // moved too far				
				histogramList.clear();			
				// @IMP! we need copy by value not copy by reference!
				// TODO: Better way to do this? 
				for (T h : histogram) 
					histogramList.add(h.copy()); 			
				Collections.sort(histogramList);	
			}
			else{ // moved by one
				RandomAccessible<T> dropSlice = Views.hyperSlice(infSrc, direction, step < 0 ? max[direction] : min[direction]);
				RandomAccessible<T> addSlice  = Views.hyperSlice(infSrc, direction, step < 0 ? min[direction] : max[direction]);				
				// drop one dimension
				long[] localMin  = new long[n - 1];
				long[] localMax  = new long[n - 1];
				
				for (int i = 0; i < n; ++i){
					if (i != direction){
						localMin[i > direction ? i - 1 : i] = min[i];
						localMax[i > direction ? i - 1 : i] = max[i];
					}

				}

				RandomAccessibleInterval<T> dropHistogram = Views.interval(dropSlice, localMin, localMax);
				RandomAccessibleInterval<T> addHistogram  = Views.interval(addSlice, localMin, localMax);

				dropElements(dropHistogram, histogramList);
				addElements(addHistogram, histogramList);

				if (debug)
					if (!isSizeCorrect(histogramList, kernel))
						return;
			}

			// get the median value 
			rDst.setPosition(cSrc);
			try{
				rDst.get().set(histogramList.get(histogramList.size()/2));
				idx++;
			}
			catch(Exception e){
				System.out.println("histogramList.size() is wrong");
				return;
			}
		}

		System.out.println("idx = " + idx);
	}
	
	/**
	 *  adds every element from histogram contained in list
	 */
	public static <T extends RealType<T> & Comparable<T>> void addElements(RandomAccessibleInterval<T> histogram, List<T> list){
		for (T h : Views.iterable(histogram)){
			int key = Collections.binarySearch(list, h);
			if (key >= 0){ // same item found
				list.add(key + 1, h.copy()); // insert after this item
			}
			else{
				key = -(key + 1);
				list.add(key, h.copy());
			}
		} 	
	}
	
	/**
	 *  removes every element from histogram contained in list
	 */
	public static <T extends RealType<T> & Comparable<T>> void dropElements(RandomAccessibleInterval<T> histogram, List<T> list){
		for (T h : Views.iterable(histogram)) {		
			int key = Collections.binarySearch(list, h);					
			try{
				list.remove(key);
			}
			catch(Exception e){
				System.out.println("Wrong value key = " + key + " specified.");
			}
		}	
	}
	
	
	/**
	 *  checks if the size of the histogram is correct
	 */
	public static <T extends RealType<T>> boolean isSizeCorrect(List<T> h, final Interval k){
		boolean res = true; 
		long totalDimensions = 1;
		for (int j = 0; j < k.numDimensions(); j++) {
			totalDimensions *= k.dimension(j);
		}

		if (h.size() > totalDimensions) { // 2D case
			System.out.println("HistogramList:" + h.size()  + " > " + totalDimensions );
			System.out.println("Something went wrong: too many elements in HistogramList");
			res = false; 
		}

		if (h.size() <  totalDimensions) {
			System.out.println("HistogramList:" + h.size()  + " > " + totalDimensions );
			System.out.println("Something went wrong: too few elements in HistogramList");
			res = false; 
		}
		return res;
	}
	
	
	/**
	 * This one checks if the cursor moved only by one pixel away
	 * @param pPos - position of the previous element
	 * @param cPos - position of the current element
	 * @param res - [0]: direction of movement (d); [1]: step
	 */
	public static void checkDist(long[] pPos, long[] cPos, long[] res){
		long n = cPos.length; // number of dimensions
		res[0] = -1;
		res[1] = 0;
		for (int d = 0; d < n; ++d){
			long dist = cPos[d] - pPos[d]; // dist > 0 if we moved forward, dist < 0 otherwise	
			if (dist != 0){ // ?moved
				if((Math.abs(dist) != 1) || (res[0] != -1)){ //?too far or ?more thatn once
					res[0] = -2;
					break;
				}
				else{
					res[0] = d; 	// set the direction of movement
					res[1] = dist;  // set the step 
				}
			}
		}
	}
	
	public static void main(String [] args){
		new ImageJ(); // to have a menu!
		
		//File file = new File("src/main/resources/Bikesgray.jpg");
		File file = new File("src/main/resources/salt-and-pepper.tif");
		// File file = new File("src/main/resources/test3D.tif");
		// File file = new File("src/main/resources/inputMedian.png");
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		
		
		final int n = img.numDimensions();
		final long[] min = new long[n];
		final long[] max = new long[n];
		
		// not super important 
		FloatType minValue = new FloatType();
		FloatType maxValue = new FloatType();
		minValue.set(0);
		maxValue.set(255);		
		Normalize.normalize(img, minValue, maxValue);
		
		ImageJFunctions.show(img);
		
		// define the size of the filter
		int zz = 1;
		// run multiple tests
		for (int jj = zz; jj <= zz; jj++) {	

			for (int d = 0; d < n; d++) {
				min[d] = -jj;
				max[d] = jj;
			}
			long inT  = System.nanoTime();
			medianFilter(img, dst, new FinalInterval(min, max));
			System.out.println("Time for jj = " + jj + " : "+ TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - inT)/1000);
			ImageJFunctions.show(dst);
		}
		
		
		//here comes filtering part 
		
		
		
		System.out.println("Doge!");
	}	
}
