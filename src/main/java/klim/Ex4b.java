package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex4b {
	
	public <T extends RealType<T> & NativeType<T>> Ex4b(){
		Img<T> img = (Img<T>) ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		Gauss.inDoubleInPlace(new double[]{1, 1}, img);
		Img <BitType> display = findAndDisplaylocalMinima(img, new ArrayImgFactory<BitType>(), new BitType());  
		ImageJFunctions.show(img);
		ImageJFunctions.show(display);		
	}
	
	public static <T extends Comparable<T>, U extends RealType<U>> Img <U> 
	findAndDisplaylocalMinima(RandomAccessibleInterval<T> src, ImgFactory<U> imageFactory, U outputType){
		Img <U> output = imageFactory.create(src, outputType);
		Interval interval = Intervals.expand(src, -1);
		src = Views.interval(src, interval);
		final Cursor<T> center = Views.iterable(src).cursor();
		final RectangleShape shape = new RectangleShape(1, true);
		for (final Neighborhood<T> localNeighborhood: shape.neighborhoods(src)){
			final T centerValue = center.next();
			boolean isMin = true; 
			for (final T value : localNeighborhood){
				if (centerValue.compareTo(value) >= 0){
					isMin = false;
					break;
				}
			}
			if (isMin){
				HyperSphere<U> hyperSphere = new HyperSphere <U> (output, center, 1);
				for (U value : hyperSphere)
					value.setOne();
			}
		}
		return output;
	}
	
	public static void main (String[] args){
		new ImageJ();
		new Ex4b();
	}
}
