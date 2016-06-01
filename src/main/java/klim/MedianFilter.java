package klim;

import java.io.File;

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
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
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
		
		while(cSrc.hasNext()){
			cSrc.localize(prev);
			cSrc.fwd();
			// cSrc.next(); // not sure if I need this one! 
			cSrc.localize(cur);
			
			// check if the cursor moved only by one step
			// movedByOne >= 0 - shows the direction of movement
			// movedByOne == -1 - initial value
			// movedByOne == -2 - moved too far
			long movedByOne = checkDist(prev, cur);
			
			if (movedByOne == -2){
				// define new boundaries of the new kernel-window
				for (int d = 0; d < n; ++d){
					min[d] = cSrc.getLongPosition(d) - kernel.dimension(d);
					max[d] = cSrc.getLongPosition(d) + kernel.dimension(d);
				}
				histogram = Views.interval(infSrc, min, max);
			
				// TODO: looks fine to search for a full median here
				// and assign the value to a new image
				
			}
			else{
				// shift by one
				// direction is given by movedByOne
				// TODO: drop elements -- add new elements 
				// search for a full median 
			}
			
		}

		
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
	public static < T extends RealType<  T > > long checkDist(long[] prev, long[] cur){
		// movedByOne >= 0 - shows the direction of movement
		// movedByOne == -1 - initial value
		// movedByOne == -2 - moved too far
		long n = cur.length;
		long movedByOne = -1;
		for (int d = 0; d < n; ++d){
			//  moved by one ?
			if (Math.abs(prev[d] - cur[d]) == 1){
				// sure we have not moved by one already?
				if (movedByOne == -1){
					// set the direction of movement
					movedByOne = d;
				}
				else{
					movedByOne = -2;
					break;
				}						
			}
			else{
				// we moved too far 
				movedByOne = -2;
				break;
			}
		}
		return movedByOne;
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
		final Img<FloatType> dst = img.copy();
		
		
		final int n = img.numDimensions();
		final long[] min = new long[n];
		final long[] max = new long[n];
		
		for (int d = 0; d < n; d++) {
			min[d] = -1;
			max[d] = 1;
		}
		
		//ImageJFunctions.show(img);
		
		//here comes filtering part 
		
		medianFilter(img, dst, new FinalInterval(min, max));
		System.out.println("Doge!");
	}
	
	
}
