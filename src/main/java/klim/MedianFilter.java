package klim;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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
		
		final long[] min = new long[n];
		final long[] max = new long[n];
		
		// histogram
		// this should be an image of the 
		
		long[] prev = new long[n]; 
		long[] cur = new long[n]; 
		IterableInterval<T> histogram;
		// the one below can be changed to LinkedList 
		// the question is what implementation is faster
		// ref: http://stackoverflow.com/questions/322715/when-to-use-linkedlist-over-arraylist?lq=1
		List<T> histogramList = new ArrayList<T>(); 
		//= new List<T>(5);
		
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
			
			long[] di = new long[2];
			
			// @TODO: Dirty! Clean this
			checkDist(prev, cur, di);
			
			direction = (int)di[0];
			step = di[1];
						
			if (direction == -2){
				// define new boundaries of the new kernel-window
				for (int d = 0; d < n; ++d){
					min[d] = cSrc.getLongPosition(d) - kernel.dimension(d);
					max[d] = cSrc.getLongPosition(d) + kernel.dimension(d);
				}
				histogram = Views.interval(infSrc, min, max);
			
				// TODO: looks fine to search for a full median here
				// and assign the value to a new image
				// looks like this part is working correctly
				histogramList.clear();
				
				for (T h : histogram) {
					histogramList.add(h);
				}
				
				Collections.sort(histogramList);
				
			}
			else{
				// shift by one
				// direction is given by movedByOne
				// TODO: drop elements -- add new elements 
				// search for a full median 
				// histogramList.clear();
				
				// so what has to be done
				// remove old values after that add new values
				
				// this view points to the part that we should cut
				// here can be the problem with the src -> outofbounds
				IterableInterval<T> oldHistogram = Views.hyperSlice(infSrc, direction, (step > 0 ? min[direction] : max[direction]));
				for (T h : oldHistogram) {
					int key = Collections.binarySearch(histogramList, h);
					histogramList.remove(key);
					
					// ?! histogramList.remove(h);
				}
				
				min[direction] = cSrc.getLongPosition(direction) - kernel.dimension(direction) + step;
				max[direction] = cSrc.getLongPosition(direction) + kernel.dimension(direction) + step;
				
				
				
				histogram = Views.interval(infSrc, min, max);
				// histogramList.clear();
				for (T h : histogram) {
					histogramList.add(h);
				}
				Collections.sort(histogramList);
				// create a view from the initial image
				// not form the histogram
				//Views.hyperSlice(histogram, direction, step > 0 ? min[direction] : max[direction]);
				
				
				
				// Collections.binarySearch(histogramList, key);
				
			}
			
			rDst.setPosition(cSrc);
			if (!histogramList.isEmpty()){
				rDst.get().set(histogramList.get(histogramList.size()/2));
				// debug 
				idx++;
			}
		}

		// System.out.println(histogramList.size());
		System.out.println("idx = " + idx);
	}
	
	// can be done using insertions 
	// instead of sorting 
	public static <T extends RealType<T>> void getMedian(final IterableInterval<T> set, final T result){
		for(T value : set){
			
		}
		
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
	
	
	public static < T extends RealType<  T > & Comparable<T> > int checkDist2(RandomAccessibleInterval<T> prev, RandomAccessibleInterval <T> cur){
		// movedByOne >= 0 - shows the direction of movement
		// movedByOne == -1 - initial value
		// movedByOne == -2 - moved too far
		// System.out.println("Moved by " + );
		int movedByOne = -1;
		
		Cursor<T> prevCursor = Views.iterable(prev).cursor();
		RandomAccess<T> curRandomAccess = cur.randomAccess();
		
		long n = cur.numDimensions();
		
		
		while(prevCursor.hasNext()){
			prevCursor.fwd();
			curRandomAccess.setPosition(prevCursor);
		}
		
//		for (int d = 0; d < n; ++d){
//			//  moved by one ?
//			if (
//					
//					Math.abs(prev. - cur[d]) == 1){
//				// sure we have not moved by one already?
//				if (movedByOne == -1){
//					// set the direction of movement
//					movedByOne = d;
//				}
//				else{
//					movedByOne = -2;
//					break;
//				}						
//			}
//			else{
//				// we moved too far 
//				movedByOne = -2;
//				break;
//			}
//		}
		return movedByOne;
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
	
	
	
	public static void main(String [] args){
		File file = new File("src/main/resources/Bikesgray.jpg");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		
		
		final int n = img.numDimensions();
		final long[] min = new long[n];
		final long[] max = new long[n];
		
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
	
	
}
