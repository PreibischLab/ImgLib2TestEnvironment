package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex2c {
 
	public <T extends Type <T>> void copy (RandomAccessible<T> src, IterableInterval<T> dst){
		Cursor<T> dstCusrsor = dst.localizingCursor();
		RandomAccess <T> srcRandomAccess = src.randomAccess();
		
		while (dstCusrsor.hasNext()){
			dstCusrsor.fwd();
			srcRandomAccess.setPosition(dstCusrsor);
			dstCusrsor.get().set(srcRandomAccess.get());
		}	
	}
	
	public Ex2c(){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		Img<FloatType> duplicate = img.factory().create(img, img.firstElement());
		copy(img, duplicate);
		ImageJFunctions.show(duplicate);
		
		RandomAccessibleInterval<FloatType> viewSrc = Views.offsetInterval(img, new long[]{100, 100}, new long[]{250, 150});
		RandomAccessibleInterval< FloatType > viewDst = Views.offsetInterval( img,new long[] { 500, 200 }, new long[]{ 250, 150 } );
		IterableInterval<FloatType> iterableDst = Views.iterable(viewDst);
		
		copy(viewSrc, iterableDst);
		
		ImageJFunctions.show(img);
		
		}
	
	public static void main (String[] args){
		new ImageJ();
		new Ex2c();
	}
}
