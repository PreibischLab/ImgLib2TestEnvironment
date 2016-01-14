package david;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MeanFilter
{
	public static <T extends NumericType<T>> void meanFilterSphere(RandomAccessibleInterval<T> src, RandomAccessibleInterval<T> dest, int radius){
		T ONE = Views.iterable(src).firstElement().createVariable();
		ONE.setOne();
		
		Cursor<T> c = Views.iterable(src).localizingCursor();
		RandomAccess<T> r = dest.randomAccess();
		while (c.hasNext()){
			c.fwd();
			
			T sum = c.get().createVariable();
			T count = c.get().createVariable();
			sum.setZero();
			count.setZero();
			
			HyperSphere<T > hs = new HyperSphere<T>(src, c, radius);
			Cursor<T > hsC = hs.localizingCursor();
			
			while (hsC.hasNext()){
				hsC.fwd();
				if (Intervals.contains(src, hsC)){
					sum.add(hsC.get());
					count.add(ONE);
					Util.printCoordinates(hsC);
				}
			}
			
			sum.div(count);
			r.setPosition(c);
			r.get().set(sum);
			
			
		}
		
	}
	
	public static <T extends NumericType<T>> void meanFilterBox(RandomAccessibleInterval<T> src, RandomAccessibleInterval<T> dest, int radius)
	{
		
		T ONE = Views.iterable(src).firstElement().createVariable();
		ONE.setOne();
		
		Cursor<T> c = Views.iterable(src).localizingCursor();
		RandomAccess<T> r = dest.randomAccess();
		
		long[] from = new long[src.numDimensions()];
		long[] to = new long[src.numDimensions()];		
		
		while (c.hasNext()){
			c.fwd();
			c.localize(from);
			c.localize(to);
			for (int i = 0; i < src.numDimensions(); i++){
				from[i] -= radius;
				to[i] += radius;
			}
			T sum = c.get().createVariable();
			T count = c.get().createVariable();
			sum.setZero();
			count.setZero();
			
			Interval tInterval = new FinalInterval(from, to);
			tInterval = Intervals.intersect(tInterval, src);
			
			
			Cursor<T> cNeighborhood = Views.interval(src, tInterval).cursor();
			
			while (cNeighborhood.hasNext()){
				cNeighborhood.fwd();
				sum.add(cNeighborhood.get());
				count.add(ONE);
				Util.printCoordinates(cNeighborhood);
			}
			
			sum.div(count);
			r.setPosition(c);
			r.get().set(sum);
						
		}
		
	}
	
	public static void main(String[] args) {
		new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		Img<FloatType> filtered = img.factory().create(img, img.firstElement());
		meanFilterBox(img, filtered, 5);
		ImageJFunctions.show(filtered);

		
		
		
	}

}
