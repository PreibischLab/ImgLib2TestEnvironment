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

public class min {
	
	public static <T extends Comparable <T>> RandomAccess<T> Globalmin (RandomAccessibleInterval<T> img){
		
		Cursor<T> c =Views.iterable(img).cursor();
		
		c.fwd();
		RandomAccess<T> m = img.randomAccess();
		m.setPosition(c);
		
		
		
		while (c.hasNext()){
			c.fwd();
			if (m.get().compareTo(c.get()) > 0){
				m.setPosition(c);
			}
		}
				
		return m;
		
		
	}
	
	

	public static void main(String[] args) {
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		
		
		long[] pos = new long[img.numDimensions()];
		RandomAccess<FloatType> m = Globalmin(img);
		m.localize(pos);
		FloatType val = m.get();
		
		
		
		System.out.println(Util.printCoordinates(pos));
		System.out.print(val);
		
		
		
		
	}
	}







