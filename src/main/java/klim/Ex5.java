package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.ExtendedRandomAccessibleInterval;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex5 {
	public Ex5(){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		RandomAccessible<FloatType> infinite1 = Views.extendValue(img, new FloatType(0));
		RandomAccessible<FloatType> infinite2 = Views.extendValue(img, new FloatType(128));
		RandomAccessible<FloatType> infinite3 = Views.extendRandom(img, 0, 255);
		RandomAccessible<FloatType> infinite4 = Views.extendMirrorSingle(img);
		RandomAccessible<FloatType> infinite5 = Views.extendMirrorDouble(img);
		RandomAccessible<FloatType> infinite6 = Views.extendPeriodic(img);
		RandomAccessible<FloatType> infinite7 = new ExtendedRandomAccessibleInterval<FloatType, Img<FloatType>>(img, 
				new OutOfBoundsConstantValueFactory<FloatType, Img<FloatType>>(new FloatType(255))); 
		
		long[] min = new long[img.numDimensions()];
		long[] max = new long[img.numDimensions()];
		
		for (int d = 0; d < img.numDimensions(); ++d){
			min[d] = -img.dimension(d) - 90;
			max[d] = img.dimension(d)*2 - 1 + 90;
		}
		
		FinalInterval interval = new FinalInterval(min, max);
		
		ImageJFunctions.show(Views.interval(infinite1, interval));
		ImageJFunctions.show(Views.interval(infinite2, interval));
		ImageJFunctions.show(Views.interval(infinite3, interval));
		ImageJFunctions.show(Views.interval(infinite4, interval));
		ImageJFunctions.show(Views.interval(infinite5, interval));
		ImageJFunctions.show(Views.interval(infinite6, interval));
		ImageJFunctions.show(Views.interval(infinite7, interval));
	}
	
	public static void main (String[] args){
		new ImageJ();
		new Ex5();
	}
}
