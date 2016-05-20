package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex6a3 {
	public Ex6a3() throws IncompatibleTypeException{
		Img <FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		double sigma = 16;
		RandomAccessible<FloatType> infiniteImg = Views.extendMirrorSingle(img);
		FinalInterval interval = Intervals.createMinMax(100, 30, 500, 250);
		RandomAccessibleInterval<FloatType> region = Views.interval(img, interval);
		Gauss3.gauss(sigma, infiniteImg, region);
		ImageJFunctions.show(img);
	}
	public static void main(String[] args) throws IncompatibleTypeException{
		new ImageJ();
		new Ex6a3();
	}
}
