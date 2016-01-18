package david;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

/**
 * 
 * For every pixel see if it is above or below a threshold
 * 
 * if above, set 1, otherwise 0 > use a binary image for this (BitType)
 *
 */
public class Thresholding
{
	
	public static <T extends Comparable<T>> void threshold(IterableInterval<T> img, RandomAccessibleInterval<BitType> result, T threshold)
	{
		Cursor<T> c = img.localizingCursor();
		RandomAccess<BitType> r = result.randomAccess();
		
		while(c.hasNext()){
			c.fwd();
			r.setPosition(c);
			if (c.get().compareTo(threshold) > 0){
				r.get().setOne();
			} else {
				r.get().setZero();
			}
		}
	}
	
	public static void main(String[] args) {
		new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		Img<BitType> thresholded = new ArrayImgFactory<BitType>().create(img, new BitType());
		threshold(img, thresholded, new FloatType(100));
		ImageJFunctions.show(thresholded);
		
		
	}

}
