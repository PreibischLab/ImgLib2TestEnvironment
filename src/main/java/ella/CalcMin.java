package ella;

import java.io.File;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Cursor;
import net.imglib2.view.Views;
import util.ImgLib2Util;
 

public class CalcMin {

	public static <T extends Comparable<T>> RandomAccess<T> findMin(RandomAccessibleInterval<T> img){
		
		final Cursor<T> c = Views.iterable(img).cursor();
		c.fwd();
		RandomAccess<T> min = img.randomAccess();
		min.setPosition(c);
		
		while ( c.hasNext()) {
			c.fwd();
			if (min.get().compareTo(c.get()) > 0){
				min.setPosition(c);
			}
		}
		
		return min;
		
	}
	
	public static void main(String[] args) {
		Img< FloatType > img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		RandomAccess<FloatType> min = findMin(img);
		
		System.out.print(Util.printCoordinates(min) + " " + min.get());
		
	}
	
}
