package david;

import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class SphereDrawing {
	
	
	public static <T extends RealType<T>> void drawSpheresCenter(RandomAccessibleInterval<T> img, long radius, double smallRadius)
	{
		Point center = new Point(img.numDimensions());
		for (int i = 0; i < img.numDimensions(); i++){
			center.setPosition(img.dimension(i) / 2, i);
		}
		
		Random rng = new Random();
		
		HyperSphere<T> hs = new HyperSphere<T>(Views.extendValue(img, Views.iterable(img).firstElement().createVariable()), center, radius);
		Cursor<T> c = hs.cursor();
		
	
		
		while (c.hasNext()){
			c.fwd();
			
			
			
			if (rng.nextDouble() < Math.pow(10, (float) - img.numDimensions())){
				T value = Views.iterable(img).firstElement().createVariable();
				value.setOne();
				value.mul(rng.nextDouble());	
				drawSphere(img, (long) (Math.ceil(rng.nextGaussian() * smallRadius)), c, value);
			}
			
		}
		
		
		
	}
	
	public static <T extends RealType<T>> void drawSphere(RandomAccessibleInterval<T> img,
																	long size,
																	Localizable center,
																	T value){

		HyperSphere<T> hs = new HyperSphere<T>(Views.extendValue(img, Views.iterable(img).firstElement().createVariable()), center, size);
		Cursor<T> c = hs.cursor();
		
		
		while (c.hasNext()){
			c.fwd();
			if (c.get().compareTo(value) < 0){
				c.get().set(value);
			}
		}
		
		
	}
	
	public static void main(String[] args) {
		Img<FloatType> img = ArrayImgs.floats(500, 500, 200);
		new ImageJ();
		drawSpheresCenter(img, 200, 6);
		ImageJFunctions.show(img);
	}

}
