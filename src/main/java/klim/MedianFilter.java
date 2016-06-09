package klim;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

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
	
	
	// TODO: write description 
	/*
	 * 
	 * */
	public static /*correct type declaration */< T extends RealType<  T > > void medianFilter(
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
		long[] prev = new long[n]; 
		long[] cur  = new long[n]; 
		
		// contains all elements of the kernel
		IterableInterval<T> histogram = Views.interval(infSrc, min, max);

		// the one below can be changed to LinkedList 
		// the question is what implementation is faster
		// ref: http://stackoverflow.com/questions/322715/when-to-use-linkedlist-over-arraylist?lq=1
		List<T> histogramList = new ArrayList<T>(); 
		
		for (T h : histogram) {
			histogramList.add(h.copy()); 
		}			
		
		Collections.sort(histogramList);	
		
		// @DEBUG: delete idx after done 
		long idx = 1;
		while(cSrc.hasNext()){
			cSrc.localize(prev);
			cSrc.fwd();
			cSrc.localize(cur);
	
			// check if the cursor moved only by one step
			// direction >= 0 - shows the direction of movement
			// direction == -1 - initial value
			// direction == -2 - moved too far
			int direction = -2;
			long step = 0; // shows the direction of movement (+1/-1)
				
			// @TODO: Dirty! Clean this
			long[] di = new long[2];
			checkDist(prev, cur, di);		
			direction = (int)di[0];
			step = di[1];
			
			// define new boundaries of the new kernel-window
			for (int d = 0; d < n; ++d){
				min[d] = cSrc.getLongPosition(d) - kernel.dimension(d)/2;
				max[d] = cSrc.getLongPosition(d) + kernel.dimension(d)/2;
			}
			
			// adjust the histogram window
			histogram = Views.interval(infSrc, min, max);
						
			if (direction == -2){ // moved too far
				
				
//				System.out.println("NEW LINE:");
//				for (int d = 0; d < n; ++d){
//					System.out.print("d[" + d + "]: ");
//					System.out.print("(" + min[d] + ", " + cSrc.getLongPosition(d) + ", " + max[d] + ") ");
//				}
//				System.out.println();
				
				
				
				// clear histogram to add new values
				histogramList.clear();			
				// @IMP! we need copy by value not copy by reference!
				// @TODO: Better way to do this? 
				for (T h : histogram) {
					histogramList.add(h.copy()); 
				}			
				// this should be the only sorting that is performed
				Collections.sort(histogramList);
		
			}
			else{ // moved by one
		
//				for (int d = 0; d < n; ++d){
//					System.out.print("d[" + d + "]: ");
//					System.out.print("(" + min[d] + ", " + cSrc.getLongPosition(d) + ", " + max[d] + ") ");
//				}
//				System.out.println();
				
				// @TODO: check that the assignment is correct
				RandomAccessible<T> dropSlice = Views.hyperSlice(infSrc, direction, step < 0 ? max[direction] : min[direction]);
				RandomAccessible<T> addSlice  = Views.hyperSlice(infSrc, direction, step < 0 ? min[direction] : max[direction]);
				
				// drop one dimension
				long[] localMin  = new long[n - 1];
				long[] localMax  = new long[n - 1];
				
				for (int i = 0; i < n; ++i){
					if (i != direction){
						localMin[i >= direction ? i - 1 : i] = min[i];
						localMax[i >= direction ? i - 1 : i] = max[i];
					}

				}

				RandomAccessibleInterval<T> dropHistogram = Views.interval(dropSlice, localMin, localMax);
				RandomAccessibleInterval<T> addHistogram  = Views.interval(addSlice, localMin, localMax);
				
//				System.out.println("histogram:");
//				for (T h : histogram) {	
//					System.out.print(h + " ");
//				}
//				System.out.println();
//				System.out.println("add      :");
//				for (T h : Views.iterable(addHistogram)) {											
//					System.out.print(h + " ");
//				}
//				System.out.println();
//				System.out.println("del      :");
//				for (T h : Views.iterable(dropHistogram)) {											
//					System.out.print(h + " ");
//				}
//				System.out.println();
//				System.out.println("list     :");
//				for (T h : histogramList) {											
//					System.out.print(h + " ");
//				}
//				System.out.println();
				
				for (T h : Views.iterable(dropHistogram)) {					
					histogramList.remove(h);
				}

				
				for (T h : Views.iterable(addHistogram)){
					histogramList.add(h.copy());
				}
				//System.out.println();		
				
				Collections.sort(histogramList);
				// System.out.println("Length after: " + histogramList.size());
			
				
//				System.out.println("(" + min[direction] + ", " + cSrc.getLongPosition(direction) + ", " + max[direction] + ")");
				// System.out.println();
				
				// histogram = Views.interval(infSrc, min, max);	
				
				
				// Collections.binarySearch(histogramList, key);
				
				if (histogramList.size() > kernel.dimension(0)*kernel.dimension(1)) {
					System.out.println("HistogramList:" + histogramList.size()  + " > " + kernel.dimension(0)*kernel.dimension(1) );
					System.out.println("Something went wrong: too many elements in HistogramList");
					
//					for (T h : histogramList) {	
//						System.out.print(h + " ");
//					}
//					System.out.println();
//					for (T h : Views.iterable(addHistogram)) {											
//						System.out.print(h + " ");
//					}
//					System.out.println();
//					for (T h : Views.iterable(dropHistogram)) {											
//						System.out.print(h + " ");
//					}
//					System.out.println();				
					
					return; 
				}
				
				if (histogramList.size() <  kernel.dimension(0)*kernel.dimension(1)) {
					System.out.println("Beep");
					return;
				}
				
			}
			
			// get the median value 
			rDst.setPosition(cSrc);
			if (!histogramList.isEmpty()){
				rDst.get().set(histogramList.get(histogramList.size()/2));
				// debug 
				idx++;
			}	
			else{
				// catch errors
				// System.out.println("Something went wrong");
			}
		}

		// System.out.println(histogramList.size());
		System.out.println("idx = " + idx);
	}
	
	/**
	 * This one checks if we are only one step away from the previous pixel
	 * @param x - previous element
	 * @param y - current element
	 * @param n - number of dimensions 	
	 */
	public static < T extends RealType<  T > > void checkDist2(long[] prev, long[] cur, long[] di){
		long n = cur.length; // number of dimensions
		di[0] = -1;
		di[1] = 0;
		for (int d = 0; d < n; ++d){
			long dist = cur[d] - prev[d]; // dist > 0 if we moved forward, dist < 0 otherwise
			//long dist = prev[d] - cur[d];
			// didn't move
			if (dist == 0){
				// do nothing
			}
			else{
				//  moved by one ? 
				if ((dist == 1) || (dist == -1)){
					// sure we have not moved by one already?
					if (di[0] == -1){
						// set the dimension of movement
						di[0] = d;
						// set the direction 
						di[1] = dist;
					}
					else{
						di[0] = -2;
						break;
					}						
				}
				else{
					// we moved too far 
					di[0] = -2;
					break;
				}
			}
		}
	}
	
	/**
	 * This one checks if we are only one step away from the previous pixel
	 * @param pPos - position of the previous element
	 * @param cPos - position of the current element
	 * @param res - [0]: directions of movement; [1]: step size
	 */
	public static < T extends RealType<  T > > void checkDist(long[] pPos, long[] cPos, long[] res){
		long n = cPos.length; // number of dimensions
		res[0] = -1;
		res[1] = 0;
		for (int d = 0; d < n; ++d){
			long dist = cPos[d] - pPos[d]; // dist > 0 if we moved forward, dist < 0 otherwise
			// didn't move
			if (dist == 0){
				// do nothing
			}
			else{
				//  moved by one ? 
				if ((dist == 1) || (dist == -1)){
					// sure we have not moved by one already?
					if (res[0] == -1){
						// set the dimension of movement
						res[0] = d;
						// set the direction 
						res[1] = dist;
					}
					else{
						res[0] = -2;
						break;
					}						
				}
				else{
					// we moved too far 
					res[0] = -2;
					break;
				}
			}
		}
	}
		
	
	
	public static void main(String [] args){
		//File file = new File("src/main/resources/Bikesgray.jpg");
		//File file = new File("src/main/resources/salt-and-pepper.tif");
		File file = new File("src/main/resources/inputMedian.png");
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
		int zz = 2;
		// run multiple tests
		for (int jj = zz; jj <= zz; jj++) {	
			for (int d = 0; d < n; d++) {
				min[d] = -jj;
				max[d] = jj;
			}
			medianFilter(img, dst, new FinalInterval(min, max));
			ImageJFunctions.show(dst);
		}
		
		
		//here comes filtering part 
		
		
		
		System.out.println("Doge!");
	}	
}
