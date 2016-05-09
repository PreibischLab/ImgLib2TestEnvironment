package david;

import java.io.File;
import java.util.PriorityQueue;

import david.DistanceSource.DistanceType;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class DistanceTransformFast {
	
	/**
	 * return a 1D view of a line along the axis specified by dimToKeep
	 * at the positions in the other dimensions specified by pos
	 * @param img
	 * @param pos position of the line in the other dimensions, the position in dimToKeep has to be left out
	 * @param dimToKeep
	 * @return
	 */
	public static <T extends RealType<T>> RandomAccessibleInterval<T> getLine(RandomAccessibleInterval<T> img, long[] pos, int dimToKeep)
	{
		RandomAccessibleInterval<T> res = img;
		int droppedDims = 0;

		for (int i = 0; i < img.numDimensions(); i++){
			if (i != dimToKeep){
				res = Views.hyperSlice(res, i-droppedDims, pos[droppedDims]);
				droppedDims++;
			}
		}
		return Views.dropSingletonDimensions(res);	
	}
	
	
	/**
	 * process a single line parallel to an axis of the image
	 * @param line
	 * @param dsType
	 */
	public static <T extends RealType<T>> void processLine(RandomAccessibleInterval<T> line, DistanceType dsType){
		
		// create a DistanceSource for each finite pixel in the current dt
		PriorityQueue<DistanceSource> distanceSources = new PriorityQueue<DistanceSource>();
		Cursor<T> c = Views.iterable(line).cursor();
		while (c.hasNext()) {
			c.fwd();
			if (c.get().getRealDouble() < c.get().getMaxValue()){
				distanceSources.add(new DistanceSource(c.getLongPosition(0), c.get().getRealDouble()));
			}
		}
		
		// process the distanceSource that will lead to an update with the smallest distance value
		// continue until there are no more updates possible
		while (distanceSources.size() > 0){
			DistanceSource ds = distanceSources.poll();
			ds.applyToImageAndGrow(line, dsType);
			if (ds.isValid()) {
				distanceSources.add(ds);
			}
		}
		
	}
	
	
	/**
	 * prepare an initial distance transform Image,
	 * setting the distance of foreground pixels to 0 and background pixels to infinity  
	 * @param mask
	 * @param dt
	 */
	public static <T extends RealType<T>> void prepareDT(RandomAccessibleInterval<BitType> mask, RandomAccessibleInterval<T> dt){
		RandomAccess<T> raDest = dt.randomAccess();
		Cursor<BitType> cSrc = Views.iterable(mask).cursor();
		while(cSrc.hasNext()){
			cSrc.fwd();
			raDest.setPosition(cSrc);
			if (cSrc.get().getInteger() == 1){
				raDest.get().setReal(0);
			} else {
				raDest.get().setReal(Double.MAX_VALUE);
			}
		}
	}
	
	/**
	 * calculate the distance transform on an image prepared by prepareDT()
	 * @param dt
	 * @param dstType
	 */
	public static <T extends RealType<T>> void calcDT(RandomAccessibleInterval<T> dt, DistanceType dstType){
		
		long [] pos = new long[dt.numDimensions()-1];
		// go through all dimensions and lines along the corresponding axis
		for (int dToKeep = 0; dToKeep < dt.numDimensions(); dToKeep ++){			
			RandomAccessibleInterval<T> hs = Views.hyperSlice(dt, dToKeep, 0);
			Cursor<T> c = Views.iterable(hs).cursor();
			while (c.hasNext()){
				c.fwd();
				c.localize(pos);
				RandomAccessibleInterval<T> line = getLine(dt, pos, dToKeep);
				processLine(line, dstType);				
				//System.out.println(Util.printCoordinates(pos));
			}
		}		
	}
	
	public static void main(String[] args) {
		new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(
				"src/main/resources/dt.png"));
		Img<BitType> thresholded = new ArrayImgFactory<BitType>().create(img,
				new BitType());
		Thresholding.threshold(img, thresholded, new FloatType(150));
		
		Img<FloatType> dt = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		//ImageJFunctions.show(getLine(img, new long[] {0, 0}, 1));
		
		prepareDT(thresholded, dt);
		//RandomAccessibleInterval<FloatType> a = getLine(img, new long[] {5}, 0);
		//RandomAccessibleInterval<FloatType> b = getLine(dt, new long[] {0}, 0);
		
		
		//ImageJFunctions.show(Views.dropSingletonDimensions(getLine(img, 1)));
		
		
		//ImageJFunctions.show(a);
		//processLine(a, DistanceType.EUCLIDEAN);
		//ImageJFunctions.show(a);
		
		
		calcDT(dt, DistanceType.EUCLIDEAN);
		
		ImageJFunctions.show(dt);
		
		/*
		RandomAccessibleInterval<FloatType> dt2 = DistanceTransform.distanceTransformSimple(thresholded);
		ImageJFunctions.show(dt2);
		*/
		
		
	}
}
