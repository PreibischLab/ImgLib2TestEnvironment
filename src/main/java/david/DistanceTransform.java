package david;

import java.io.File;
import david.Thresholding;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;

/**
 * Based on a binary image (from thresholding), compute for every bright pixel the closest distance to a 0-pixel 
 * 
 *
 */
public class DistanceTransform
{
	
	public static RandomAccessibleInterval<FloatType> distanceTransformSimple(RandomAccessibleInterval<BitType> src){
		RandomAccessibleInterval<FloatType> result = new ArrayImgFactory<FloatType>().create(src, new FloatType());
		
		Cursor<BitType> c = Views.iterable(src).localizingCursor();
		RandomAccess<FloatType> ra = result.randomAccess();
		
		while (c.hasNext()){
			c.fwd();
			
			Cursor<BitType> c2 = Views.iterable(src).localizingCursor();
			
			double dist = Double.MAX_VALUE;
			while (c2.hasNext()){
				c2.fwd();
				
				if (c2.get().getInteger() == 1 && Util.distance(c, c2) < dist){
					dist = Util.distance(c, c2);
					
				}
			}
						
			ra.setPosition(c);
			ra.get().set((float)dist);			
			
		}
		
		
		
		return result;
	}
	
	
	public static void setLocation(RandomAccessibleInterval<ShortType> locations, Localizable locationtoSet, Localizable pixel){
		RandomAccess<ShortType> ra = locations.randomAccess();
		long[] pos = new long[pixel.numDimensions() + 1];
		for ( int i = 0; i < pixel.numDimensions(); i++){
			pos[i] = pixel.getLongPosition(i);
		}
		for (int i = 0; i < locationtoSet.numDimensions(); i++){
			pos[pixel.numDimensions()] = i;
			ra.setPosition(pos);
			ra.get().set((short)locationtoSet.getIntPosition(i));
		}
		
		
	}
	
	public static Point getLocation(RandomAccessibleInterval<ShortType> locations, Localizable pixel){
		RandomAccess<ShortType> ra = locations.randomAccess();
		long[] pos = new long[pixel.numDimensions() + 1];
		for ( int i = 0; i < pixel.numDimensions(); i++){
			pos[i] = pixel.getLongPosition(i);
		}
		long[] res = new long[(int) locations.dimension(pixel.numDimensions())];
		for(int i = 0; i < res.length; i++){
			pos[pixel.numDimensions()] = i;
			ra.setPosition(pos);
			res[i] = ra.get().get();
		}
		return new Point(res);
		
	}
	
	public static RandomAccessibleInterval<FloatType> distanceTransform(RandomAccessibleInterval<BitType> src){
		
		RandomAccessibleInterval<FloatType> result = new ArrayImgFactory<FloatType>().create(src, new FloatType());
		RandomAccessibleInterval<FloatType> resultT1 = new ArrayImgFactory<FloatType>().create(src, new FloatType());
		
		long[] dimensions = new long[src.numDimensions()+1];
		for (int i = 0; i < src.numDimensions(); i++){
			dimensions[i] = src.dimension(i);
		}
		dimensions[src.numDimensions()] = src.numDimensions();
				
		RandomAccessibleInterval<ShortType> locations = new ArrayImgFactory<ShortType>().create(dimensions, new ShortType((short)-1));
		RandomAccessibleInterval<ShortType> locationsT1 = new ArrayImgFactory<ShortType>().create(dimensions, new ShortType((short)-1));
				
		Cursor<BitType> c = Views.iterable(src).localizingCursor();
		RandomAccess<FloatType> ra = result.randomAccess();
		
		long[] myMinusOne = new long[src.numDimensions()];
		for (int i = 0; i < src.numDimensions(); i++){myMinusOne[i] = -1;}
		Point minusOnePoint = new Point(myMinusOne);
		
		
		while(c.hasNext()){
			c.fwd();
			ra.setPosition(c);
			if (c.get().getInteger() == 1){
				ra.get().set(0);
				setLocation(locations, ra, ra);
			} else {
				ra.get().set(Float.MAX_VALUE);
				
				setLocation(locations, minusOnePoint, ra);
			}
		}
		

		boolean stopped = false;
		
//		long[] max = new long[result.numDimensions()];
//		long[] min = new long[result.numDimensions()];
//		result.max(max);
//		result.min(min);
//		double max_distance = Util.distance(min, max);
//		
//		System.out.println(max_distance);
		
		ImageStack stack = new ImageStack( (int)src.dimension(0), (int)src.dimension(1));
		
		stack.addSlice( ImageJFunctions.wrap( result, "").getProcessor() );
		
		int i =0;
		
		while ( !stopped){			
			
			stopped = dtStep(result, resultT1, locations, locationsT1);
			
//			final Cursor<FloatType> c1 = Views.iterable(result).cursor();
//			final Cursor<FloatType> c2 = Views.iterable(resultT1).cursor();
//			
//			boolean change = false;
//			while (c1.hasNext()){
//				if (c1.next().get() != c2.next().get()){
//					change = true;
//				}
//			}
			
			DavidTest.copyImage(resultT1, Views.iterable(result));
			DavidTest.copyImage(locationsT1, Views.iterable(locations));
			
			stack.addSlice( ImageJFunctions.wrap( result, "").getProcessor() );
			
//			System.out.println( "iteration " + (++i) + "stopped:" + stopped + " change:" + change );
			System.out.println( "iteration " + (++i));
			
		}

		new ImagePlus( "illustration", stack ).show();
		return result;
		
		
	}
	
	public static boolean dtStep(RandomAccessibleInterval<FloatType> src, RandomAccessibleInterval<FloatType> result,
								RandomAccessibleInterval<ShortType> locations, RandomAccessibleInterval<ShortType> locations2){
		Cursor<FloatType> cSrc = Views.iterable(src).localizingCursor();
		RandomAccess<FloatType> raDest = result.randomAccess();
		
		long[] myMinusOne = new long[src.numDimensions()];
		for (int i = 0; i < src.numDimensions(); i++){myMinusOne[i] = -1;}
		Point minusOnePoint = new Point(myMinusOne);
		
		boolean stopped = true;
		boolean changed = false;
		
		while(cSrc.hasNext()){
			cSrc.fwd();
			raDest.setPosition(cSrc);
			
			RandomAccess<FloatType> minRA = getMinDistanceIn4NeighborhoodRA(cSrc, src, locations);
			
			for (int i = 0; i < cSrc.numDimensions(); i++){
				if (cSrc.getIntPosition(i) != minRA.getIntPosition(i)) {changed = true;}
			}
			
			
			if ( getLocation(locations, minRA).getIntPosition(0) != -1 ){
			
				setLocation(locations2, getLocation(locations, minRA), cSrc);
			
				raDest.get().set((float)Util.distance(cSrc, getLocation(locations, minRA)));
			} else {
				setLocation(locations2, minusOnePoint, cSrc);
				raDest.get().set(Float.MAX_VALUE);
			}
			
			
			if (cSrc.get().get() == Float.MAX_VALUE){
				
				stopped = false;
			}
		}
		
		return stopped;
		//return !changed;
	}
	
	public static FloatType getMinDistanceInNeighborhood(Localizable x, RandomAccessibleInterval<FloatType> src){
		long[] from = new long[src.numDimensions()];
		long[] to = new long[src.numDimensions()];
		x.localize(from);
		x.localize(to);
		for (int i = 0; i < src.numDimensions(); i++){
			from[i]--;
			to[i]++;
		}
		Interval tInterval = new FinalInterval(from, to);
		
		
		RandomAccess<FloatType> minRA = MinSearch.findMinWithSmallestDistanceToCenter(Views.interval(Views.extendValue(src, new FloatType(Float.MAX_VALUE)), tInterval), x);
		FloatType min = new FloatType((float)Util.distance(x, minRA));
		min.add(minRA.get());
		return min;
		
	}
	
	public static RandomAccess<FloatType> getMinDistanceIn4NeighborhoodRA(Localizable x, RandomAccessibleInterval<FloatType> src,
			RandomAccessibleInterval<ShortType> locations){
		RandomAccess<FloatType> minRA = Views.extendValue(src, new FloatType(Float.MAX_VALUE)).randomAccess();
		minRA.setPosition(x);
		
		RandomAccess<FloatType> r = Views.extendValue(src, new FloatType(Float.MAX_VALUE)).randomAccess();
		
		long[] pos = new long[src.numDimensions()];
		
		for (int i = 0; i < src.numDimensions(); i++){
			x.localize(pos);
			if (pos[i] < src.dimension(i)) {
				pos[i]++;
				r.setPosition(pos);
				if (r.get().get() < minRA.get().get()) {
					minRA.setPosition(r);
				} else if ((r.get().get() == minRA.get().get())&& r.get().get() != Float.MAX_VALUE) {
					if (getLocation(locations, r).getIntPosition(0) != -1 && Util.distance(getLocation(locations, r),
							x) < Util.distance(getLocation(locations, minRA), x)) {
						minRA.setPosition(r);

					}
				}

			}
	

			x.localize(pos);
			if (pos[i] > 0){
				pos[i]--;
				r.setPosition(pos);
				if (r.get().get() < minRA.get().get()) {
					minRA.setPosition(r);
				} else if ((r.get().get() == minRA.get().get())&& r.get().get() != Float.MAX_VALUE) {
					if (getLocation(locations, r).getIntPosition(0) != -1 && Util.distance(getLocation(locations, r),
							x) < Util.distance(getLocation(locations, minRA), x)) {
						minRA.setPosition(r);

					}
				}
			}
			
		}
		
		return minRA;
		
		
		
	}
	
	
	public static RandomAccess<FloatType> getMinDistanceInNeighborhoodRA(Localizable x, RandomAccessibleInterval<FloatType> src){
		long[] from = new long[src.numDimensions()];
		long[] to = new long[src.numDimensions()];
		x.localize(from);
		x.localize(to);
		for (int i = 0; i < src.numDimensions(); i++){
			from[i]--;
			to[i]++;
		}
		Interval tInterval = new FinalInterval(from, to);
		tInterval = Intervals.intersect(tInterval, src);
		
		
		RandomAccess<FloatType> minRA = MinSearch.findMinWithSmallestDistanceToCenter(Views.interval(Views.extendValue(src, new FloatType(Float.MAX_VALUE)), tInterval), x);
		return minRA;
		
	}
	
	
	public static void main(String[] args) {
		new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		Img<BitType> thresholded = new ArrayImgFactory<BitType>().create(img, new BitType());
		Thresholding.threshold(img, thresholded, new FloatType(200));
		ImageJFunctions.show(thresholded);
//		RandomAccessibleInterval<FloatType> res = distanceTransformSimple(thresholded);
//		ImageJFunctions.show(res);
		RandomAccessibleInterval<FloatType> res2 = distanceTransform(thresholded);
		ImageJFunctions.show(res2);
		
		
	}

}
