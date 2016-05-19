package klim;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import util.ImgLib2Util;

public class Ex3a3 {

	
	public <T extends RealType<T> & NativeType <T>> Ex3a3(){
		Img<T> img = (Img<T>) ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		Point locMin = new Point(img.numDimensions());
		Point locMax = new Point(img.numDimensions());
		computeMinMaxLocation(img, locMin, locMax);
		
		System.out.println("loc min: " + locMin);
		System.out.println("loc max: " + locMax);
	}
	
	public <T extends Comparable<T> & Type<T>> void computeMinMaxLocation(final IterableInterval<T> input, final Point minLoc, final Point maxLoc){
		final Cursor<T> cursor = input.cursor();
		T type = cursor.next();
		T min = type.copy();
		T max = type.copy();
		
		while(cursor.hasNext()){
			type = cursor.next();
			if (type.compareTo(max) > 0){
				max.set(type);
				maxLoc.setPosition(cursor);
			}
			if (type.compareTo(min) < 0){
				min.set(type);
				minLoc.setPosition(cursor);
			}	
		}
		
		
	} 
	
	public static void main(String[] args){
		new Ex3a3();
	}
}
