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
	private static final boolean typeList = false;
	
	public static < T extends RealType<  T > & Comparable<T> > void medianFilter(
			final RandomAccessibleInterval< T > src, final RandomAccessibleInterval< T > dst, final int[] kernelDim )
	{
		medianFilter(Views.extendMirrorSingle(src), src, dst, kernelDim);
	}
	
	// TODO: write description 
	/*
	 * 
	 * */
	// TODO:  /*correct type declaration */
	public static < T extends RealType<  T > & Comparable<T> > void medianFilter(
			final RandomAccessible< T > infSrc, final Interval srcInterval, final RandomAccessibleInterval< T > dst, final int[] kernelDim){

		final RandomAccessibleInterval<T> src = Views.interval(infSrc, srcInterval);
		final Cursor<T> cSrc = Views.iterable(src).localizingCursor();
		final RandomAccess<T> rDst = dst.randomAccess();
		
		final int n = src.numDimensions();
		// store kernel boundaries
		final long[] min = new long[n];
		final long[] max = new long[n];
		
		final int[] kernelHalfDim = new int[n]; 
		for (int d = 0; d < n; ++d){
			if ( kernelDim[d]%2 == 0 )
				throw new RuntimeException("kernel not odd");
			kernelHalfDim[d] = kernelDim[d]/2; // store dim/2
		}

		// store previous/current position of cursor
		final long[] pPos = new long[n]; 
		final long[] cPos = new long[n]; 
		cSrc.localize(cPos);
		updateKernelMinMax(min, max, cPos, kernelHalfDim, n);
		
		// contains all elements of the kernel
		IterableInterval<T> histogram = Views.interval(infSrc, min, max);
		// List<T> histogramList = new ArrayList<T>(); 
		// List<T> histogramList = new LinkedList<T>(); 
		
		List<T> histogramList  = typeList ? new LinkedList<T>() : new ArrayList<T>(); 	
		addAll(histogram, histogramList);

		final long[] localMin  = new long[n - 1];
		final long[] localMax  = new long[n - 1];

		// @DEBUG: delete idx after done 
		long idx = 1;
		while(cSrc.hasNext()){
			cSrc.localize(pPos);
			cSrc.fwd();
			cSrc.localize(cPos);

			final long checkDist = checkDist(pPos, cPos, n);	 // check if the cursor moved only by one step
			 // System.out.println("x: " + min[0] + " "  + cSrc.getLongPosition(0)  + " " + max[0] + " " + "y: " + min[1] + " " + cSrc.getLongPosition(1) + " " + max[1] );
						
			if (checkDist == 0){ // moved too far
				// define new boundaries of the new kernel-window
				updateKernelMinMax(min, max, cPos, kernelHalfDim, n);
				// adjust the histogram window
				histogram = Views.interval(infSrc, min, max);
				histogramList.clear();					
				addAll(histogram, histogramList);
			}
			else{ // moved by one
				final int dim = (int)Math.abs(checkDist) - 1; 	 //
				final long step = (long) Math.signum(checkDist); // shows the direction of movement (+1/-1)
	
				min[ dim ] += step;
				max[ dim ] += step;
				
				final RandomAccessible<T> dropSlice = Views.hyperSlice(infSrc, dim, step < 0 ? max[dim] + 1: min[dim] - 1);
				final RandomAccessible<T> addSlice  = Views.hyperSlice(infSrc, dim, step < 0 ? min[dim] : max[dim]);		
				
				// drop one dimension
				for (int i = 0; i < n; ++i){
					if (i != dim){
						localMin[i > dim ? i - 1 : i] = min[i];
						localMax[i > dim ? i - 1 : i] = max[i];
					}

				}
				
				final RandomAccessibleInterval<T> dropHistogram = Views.interval(dropSlice, localMin, localMax);
				final RandomAccessibleInterval<T> addHistogram  = Views.interval(addSlice, localMin, localMax);
				
				if (typeList){
					dropElements(dropHistogram, histogramList);
					addElements(addHistogram, histogramList);
				}
				else{
					dropElements2(dropHistogram, histogramList);
					addElements2(addHistogram, histogramList);
				}

				if (debug)
					if (!isSizeCorrect(histogramList, kernelDim))
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
	 * copy elements from the histogram to list 
	 * */
	public final static <T extends RealType<T> & Comparable<T>> void addAll(final IterableInterval<T> histogram, final List<T> list){
		for (final T h : histogram) 
			list.add(h.copy()); 			
		Collections.sort(list);	
	}
	
	public final static void updateKernelMinMax(final long[] min, final long[] max, final long[] position, final int[] kernelHalfDim, int n){
		for (int d = 0; d < n; ++d){
			min[d] = position[d] - kernelHalfDim[d];
			max[d] = position[d] + kernelHalfDim[d];
		}
	}
	
	/**
	 *  adds every element from histogram contained in list
	 */
	public final static <T extends RealType<T> & Comparable<T>> void addElements(RandomAccessibleInterval<T> histogram, List<T> list){
		// sort histogram first
		// then iterate over both at the same time and insert into list if necessary or increase index in sorted histogram or list
		// if concurrentModification exception just make a new list, ArrayList would work as you only add at the end
		for (T h : Views.iterable(histogram)){
			int key = Collections.binarySearch(list, h); // LinkedList is not RandomAccess, that's why this is slow
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
	 *  2nd implementation
	 *  adds every element from histogram contained in list
	 */
	public final static <T extends RealType<T> & Comparable<T>> void addElements2(RandomAccessibleInterval<T> histogram, List<T> list){
		
		List<T> histogramArray = new ArrayList<T>(); 
		List<T> returnArray = new ArrayList<T>(); 

		addAll(Views.iterable(histogram), histogramArray);

		int i = 0; // histogramArray index
		int j = 0; // list index

		int histogramArraySize = histogramArray.size();
		int listSize = list.size();

		while((i != histogramArraySize)  && (j != listSize)){
			if (histogramArray.get(i).compareTo(list.get(j)) > 0){
				returnArray.add(list.get(j++));
			} 
			else{
				if (histogramArray.get(i).compareTo(list.get(j)) < 0){
					returnArray.add(histogramArray.get(i++));
				}
				else{ // equality case
					returnArray.add(list.get(j++));
					returnArray.add(histogramArray.get(i++));
				}
			}

			// flush the rest
			if (i == histogramArraySize)
				while(j < listSize)
					returnArray.add(list.get(j++));
			if (j == listSize)
				while(i < histogramArraySize)
					returnArray.add(histogramArray.get(i++));		
		}

		list.clear();	
		for(final T h : returnArray) 
			list.add(h.copy());
	}
	
	/**
	 *  removes every element from histogram contained in list
	 */
	public final static <T extends RealType<T> & Comparable<T>> void dropElements(final RandomAccessibleInterval<T> histogram, final List<T> list){
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
	
	public final static <T extends RealType<T> & Comparable<T>> void dropElements2(final RandomAccessibleInterval<T> histogram, final List<T> list){
		
		List<T> histogramArray = new ArrayList<T>(); 
		addAll(Views.iterable(histogram), histogramArray);
	
		int i = 0; // histogramArray index
		int j = 0; // list index
		
		int histogramArraySize = histogramArray.size();
		while(i != histogramArraySize){ // iterate while we do not go through all elements in histogram
			if (histogramArray.get(i).compareTo(list.get(j)) == 0){
				list.remove(j);
				i++; 
			}
			else
				j++; 
		}			
	}
	
	
	
	/**
	 *  checks if the size of the histogram is correct
	 */
	public static <T extends RealType<T>> boolean isSizeCorrect(List<T> h, final int[] k){
		boolean res = true; 
		long totalDimensions = 1;
		for (int j = 0; j < k.length; j++) {
			totalDimensions *= k[j];
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
	 * @param n - num of dimensions
	 * @return d+1 in which it moved, or zero if jumped to far
	 */
	public static long checkDist(final long[] pPos, final long[] cPos, final int n){
		long dim = -1;
		long dir = 0;

		for (int d = 0; d < n; ++d){
			final long dist = cPos[d] - pPos[d]; // dist > 0 if we moved forward, dist < 0 otherwise
			
			if (dist != 0){ // ?moved
				if((Math.abs(dist) != 1) || (dim != -1)){ //?too far or ?more thatn once
					return 0;
				}
				else{
					dim = d; 	// set the direction of movement
					dir = dist;  // set the step 
				}
			}
		}

		return (dim + 1)*dir;
	}
	
	public static void main(String [] args){
		new ImageJ(); // to have a menu!
		
		// File file = new File("src/main/resources/Bikesgray.jpg");
		File file = new File("src/main/resources/salt-and-pepper.tif");
		//File file = new File("src/main/resources/noisyWoman.png");
		// File file = new File("src/main/resources/test3D.tif");
		// File file = new File("src/main/resources/inputMedian.png");
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
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
		
		//ImageJFunctions.show(img);
		
		// define the size of the filter
		int zz = 3;
		// run multiple tests
		for (int jj = zz; jj <= zz; jj += 2) {	

			for (int d = 0; d < n; d++) {
				min[d] = -jj;
				max[d] = jj;
			}
			long inT  = System.nanoTime();
			//final RandomAccessible< T > infSrc, final Interval srcInterval, final RandomAccessibleInterval< T > dst, final int[] kernelDim);
			//medianFilter(img, dst, new FinalInterval(min, max));
			medianFilter(img, dst, new int[]{jj, jj});
			System.out.println("kernel = " + jj + "x" + jj + " : "+ TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - inT)/1000.0);
			ImageJFunctions.show(dst);
		}
		
		
		//here comes filtering part 
		
		
		
		// System.out.println("Doge!");
	}	
}
