package varun;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.RealSum;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MinMax {

public static <T extends Comparable<T>> RandomAccess<T> LocalMax(RandomAccessibleInterval <T> img){

long startposx=0;
long startposy=0;
final long CellSizex=img.numDimensions()/2;
final long CellSizey=img.numDimensions()/2;
RandomAccess<T> localm = max(img);
while(startposx<=img.numDimensions()-CellSizex && startposy<=img.numDimensions()-CellSizey)

{

RandomAccessibleInterval< T > viewSource = Views.offsetInterval( img, new long[] { startposx, startposy }, new long[]{ startposx+CellSizex, startposy+CellSizey } );

localm = max(viewSource);

startposx+=CellSizex;

startposy+=CellSizey;

return localm;


}

return localm;

}
	
	public static <T extends Comparable<T>> RandomAccess<T> max(RandomAccessibleInterval<T> img){
		
		Cursor<T> c = Views.iterable(img).cursor();
		c.fwd();
		RandomAccess<T> m = img.randomAccess();
		m.setPosition(c);
		
		
		
		while (c.hasNext()){
			c.fwd();
			if (m.get().compareTo(c.get()) < 0){
				m.setPosition(c);
			}
		}
				
		return m;
		
	}
	
	
	


	public static void main(String[] args) {
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		
		long[] pos = new long[img.numDimensions()];
		RandomAccess<FloatType> m = max(img);
		m.localize(pos);
		FloatType val = m.get().copy();
		
		System.out.println(Util.printCoordinates(pos));
		System.out.print(val);

RandomAccess<FloatType> localm=LocalMax(img);
localm.localize(pos);
		FloatType localval = localm.get().copy();
		
		System.out.println(Util.printCoordinates(pos));
		System.out.print(localval);
		
		System.out.print("Hi"+img.numDimensions());

	}

}
