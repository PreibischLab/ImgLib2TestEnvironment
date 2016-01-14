package david;

import java.io.File;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MinSearch {

	public static <T extends Comparable<T>> RandomAccess<T> findMin(RandomAccessibleInterval<T> img) {

		Cursor<T> c = Views.iterable(img).cursor();
		RandomAccess<T> min = img.randomAccess();
		c.fwd();
		min.setPosition(c);

		while (c.hasNext()) {
			c.fwd();
			if (c.get().compareTo(min.get()) < 0) {
				min.setPosition(c);
			}
		}

		return min;
	}
	
	public static <T extends FloatType> RandomAccess<T> findMinWithSmallestDistanceToCenter(RandomAccessibleInterval<T> img, Localizable center) {
		Cursor<T> c = Views.iterable(img).cursor();
		RandomAccess<T> min = img.randomAccess();
		c.fwd();
		min.setPosition(c);

		while (c.hasNext()) {
			c.fwd();
			if (c.get().getRealDouble() + Util.distance(c, center)  < min.get().getRealDouble() + Util.distance(min, center)  ) {
				min.setPosition(c);
			}
		}

		return min;
	}
	
	public static void main(String[] args) {
		
		//new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		
		long[] pos = new long[img.numDimensions()];
		RandomAccess<FloatType> m = findMin(img);
		m.localize(pos);
		FloatType val = m.get().copy();
		
		System.out.println(Util.printCoordinates(pos));
		System.out.print(val);

	}

}
