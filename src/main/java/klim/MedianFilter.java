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
//		@DEBUG:		
//		long[] x = new long[2];
//		kernel.dimensions(x);
//		
//		System.out.println(x[0] + " " + x[1]);

		final RandomAccessible<T> infSrc = Views.extendMirrorSingle(src);
		final Cursor<T> cSrc = Views.iterable(src).localizingCursor();
		final RandomAccess<T> rDst = dst.randomAccess();
		
		final int n = src.numDimensions();
		// store kernel boundaries
		final long[] min = new long[n];
		final long[] max = new long[n];
		// store previous/current position of cursor
		long[] prev = new long[n]; 
		long[] cur  = new long[n]; 
		
		// contains all elements of the kernel
		IterableInterval<T> histogram = Views.interval(infSrc, min, max);;
		
		// the one below can be changed to LinkedList 
		// the question is what implementation is faster
		// ref: http://stackoverflow.com/questions/322715/when-to-use-linkedlist-over-arraylist?lq=1
		List<T> histogramList = new ArrayList<T>(); 
		
		// PRE-PROCESSING
		// perform first step to set all variables
		if (cSrc.hasNext()){
			cSrc.fwd();		
			for (int d = 0; d < n; ++d){
				min[d] = cSrc.getLongPosition(d) - kernel.dimension(d)/2;
				max[d] = cSrc.getLongPosition(d) + kernel.dimension(d)/2;
			}
			
			histogram = Views.interval(infSrc, min, max);
			for (T h : histogram) {
				histogramList.add(h.copy()); 
			}			
			Collections.sort(histogramList);			
		}
		
		
		// @DEBUG:
		long idx = 0;
		
		while(cSrc.hasNext()){
			// ?? should not the first value through an error? 
			cSrc.localize(prev);
			cSrc.fwd();
			// cSrc.next(); // not sure if I need this one! 
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
						
			if (direction == -2){ // moved too far
				// define new boundaries of the new kernel-window
				// @DEBUG: Passed!
				//System.out.println(prev[0] + " " + cur[0] + " " + prev[1] + " " + cur[1]);
				for (int d = 0; d < n; ++d){
					min[d] = cSrc.getLongPosition(d) - kernel.dimension(d)/2;
					max[d] = cSrc.getLongPosition(d) + kernel.dimension(d)/2;
				}
				
				// @DEBUG: Passed!
				// System.out.println(cur[0] + " " + cur[1]);
				// System.out.println(min[0] + " " + max[0] + " " + min[1] + " " + max[1]);
				
				histogram = Views.interval(infSrc, min, max);
				// Clear histogram to add new values
				histogramList.clear();
				// @DEBUG: Passed!
				// System.out.println(histogramList.size());				
				// @IMP! we need copy by value not copy by reference!
				// @TODO: Better way to do this? 
				for (T h : histogram) {
					histogramList.add(h.copy()); 
				}
				
				// System.out.println("Length before: " + histogramList.size());
				// @DEBUG: Passed!			
//				if (idx % 3000 == 0){
//					for (int i = 0; i < histogramList.size(); i++) {
//						System.out.print(histogramList.get(i) + " ");
//					}
//					System.out.println();
//				}
				
				Collections.sort(histogramList);
//				// @DEBUG: Passsed!				
//				if (idx % 3000 == 0){
//					for (int i = 0; i < histogramList.size(); i++) {
//						System.out.print(histogramList.get(i) + " ");
//					}
//					System.out.println();
//				}
				
			}
			else{ // moved by one

				RandomAccessible<T> dropSlice = Views.hyperSlice(infSrc, direction, step < 0? min[direction] : max[direction]);
				RandomAccessible<T> addSlice  = Views.hyperSlice(infSrc, direction, step > 0? min[direction] : max[direction]);
				
				// @DEBUG: Passed
//				System.out.println("DropD = " + dropSlice.numDimensions());
//				System.out.println("AddD  = " +  addSlice.numDimensions());
				
				// drop one dimension
				long[] localMin  = new long[n - 1];
				long[] localMax  = new long[n - 1];
				
				// ??? is this one working properly
				// for 2D : yes
				for (int i = 0; i < n; ++i){
					if (i == direction){
						continue;
					}
					localMin[i >= direction ? i - 1 : i] = min[i];
					localMax[i >= direction ? i - 1 : i] = max[i];
				}
				
				// @DEBUG: Passed
// 				System.out.println(localMin[0] + " " + localMax[0]);

				RandomAccessibleInterval<T> dropHistogram = Views.interval(dropSlice, localMin, localMax);
				RandomAccessibleInterval<T> addHistogram  = Views.interval(addSlice, localMin, localMax);
				
//				// @DEBUG: Passed!
//				if (idx % 3000 == 0){
//					for (T h : Views.iterable(dropHistogram)) {											
//						System.out.print(h + " ");
//					}
//					System.out.println();
//				}
				
//				// @Debug
//				System.out.println("DropD = " + dropHistogram.numDimensions());
//				System.out.println("AddD = " +  addHistogram.numDimensions());
//				System.out.println("DropD = " + dropHistogram.dimension(0));
//				System.out.println("AddD = " +  addHistogram.dimension(0));				
				
				// @DEBUG: 
//				if (idx % 3000 == 0){
//					for (T h : histogramList) {											
//						System.out.print(h + " ");
//					}
//					System.out.println();
//				}		
				
				// System.out.println(histogramList.size());

				System.out.println("Histogram     : ");
				for (T h : histogram) {											
					System.out.print(h + " ");
				}
				System.out.println();
				
				System.out.println("Histogram List: ");
				for (T h : histogramList) {											
					System.out.print(h + " ");
				}
				System.out.println();
				
				for (T h : Views.iterable(dropHistogram)) {					
					// System.out.println(h);
					
//					int key = Collections.binarySearch(histogramList, h);
//					System.out.println("key = " + key);
//					histogramList.remove(key);
					
					//histogramList.remove(h);
					System.out.print(h + " ");
					histogramList.remove(h);
				}
				
				System.out.println();
				for (T h : histogramList) {	
					System.out.print(h + " ");
				}
				System.out.println();				
				
				// System.out.println("Length before: " + histogramList.size());
				
				for (T h : Views.iterable(addHistogram)){
					System.out.print(h + " ");
					histogramList.add(h.copy());
				}
				System.out.println();		
				
				Collections.sort(histogramList);
				// System.out.println("Length after: " + histogramList.size());
				
//				@TODO: Keep this one!				
				min[direction] = cSrc.getLongPosition(direction) - kernel.dimension(direction)/2 + step;
				max[direction] = cSrc.getLongPosition(direction) + kernel.dimension(direction)/2 + step;
				
				histogram = Views.interval(infSrc, min, max);
				
//				min[direction] = cSrc.getLongPosition(direction)  + step;
//				max[direction] = cSrc.getLongPosition(direction)  + step;
//				
				//System.out.println(min[direction] + " " + max[direction]);
				
				
				
				// Collections.binarySearch(histogramList, key);
				
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
				 System.out.println("Something went wrong");
			}
			
			if (histogramList.size() > 9) return; 
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
	public static < T extends RealType<  T > > void checkDist(long[] prev, long[] cur, long[] di){
		// movedByOne >= 0 - shows the direction of movement
		// movedByOne == -1 - initial value
		// movedByOne == -2 - moved too far
		long n = cur.length;
		di[0] = -1;
		for (int d = 0; d < n; ++d){
			long dist = prev[d] - cur[d];
			
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
		
	
	
	public static void main(String [] args){
		//File file = new File("src/main/resources/Bikesgray.jpg");
		File file = new File("src/main/resources/salt-and-pepper.tif");
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
		
		for (int d = 0; d < n; d++) {
			min[d] = -1;
			max[d] = 1;
		}
		
		ImageJFunctions.show(img);
		
		//here comes filtering part 
		
		medianFilter(img, dst, new FinalInterval(min, max));
		ImageJFunctions.show(dst);
		
		System.out.println("Doge!");
	}

	// can be done using insertions 
	// instead of sorting 
	public static <T extends RealType<T>> void getMedian(final IterableInterval<T> set, final T result){
		for(T value : set){
			
		}
		
	}
	
//	public static <T extends Comparable<T>, U extends RealType<U>> Img <U> 
//	findAndDisplaylocalMaxima(RandomAccessibleInterval<T> src, ImgFactory<U> imageFactory, U outputType){
//		Img <U> output = imageFactory.create(src, outputType);
//		Interval interval = Intervals.expand(src, -1);
//		src = Views.interval(src, interval);
//		final Cursor<T> center = Views.iterable(src).cursor();
//		final RectangleShape shape = new RectangleShape(1, true);
//		for (final Neighborhood<T> localNeighborhood: shape.neighborhoods(src)){
//			final T centerValue = center.next();
//			boolean isMax = true; 
//			for (final T value : localNeighborhood){
//				if (centerValue.compareTo(value) <= 0){
//					isMax = false;
//					break;
//				}
//			}
//			if (isMax){
//				HyperSphere<U> hyperSphere = new HyperSphere <U> (output, center, 1);
//				for (U value : hyperSphere)
//					value.setOne();
//			}
//		}
//		return output;
//	}
	
//	public static < T extends RealType<  T > & Comparable<T> > int checkDist2(RandomAccessibleInterval<T> prev, RandomAccessibleInterval <T> cur){
//	// movedByOne >= 0 - shows the direction of movement
//	// movedByOne == -1 - initial value
//	// movedByOne == -2 - moved too far
//	// System.out.println("Moved by " + );
//	int movedByOne = -1;
//	
//	Cursor<T> prevCursor = Views.iterable(prev).cursor();
//	RandomAccess<T> curRandomAccess = cur.randomAccess();
//	
//	long n = cur.numDimensions();
//	
//	
//	while(prevCursor.hasNext()){
//		prevCursor.fwd();
//		curRandomAccess.setPosition(prevCursor);
//	}
//	
////	for (int d = 0; d < n; ++d){
////		//  moved by one ?
////		if (
////				
////				Math.abs(prev. - cur[d]) == 1){
////			// sure we have not moved by one already?
////			if (movedByOne == -1){
////				// set the direction of movement
////				movedByOne = d;
////			}
////			else{
////				movedByOne = -2;
////				break;
////			}						
////		}
////		else{
////			// we moved too far 
////			movedByOne = -2;
////			break;
////		}
////	}
//	return movedByOne;
//}
	
}
