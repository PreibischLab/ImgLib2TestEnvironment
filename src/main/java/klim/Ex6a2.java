package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex6a2 {

	public Ex6a2() throws IncompatibleTypeException{
		Img <FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		double[] sigma = new double[img.numDimensions()];
		for (int d = 0; d < img.numDimensions(); ++d)
			sigma[d] = 8;
		
		RandomAccessible<FloatType> infiniteImg = Views.extendValue(img, new FloatType());
		Gauss3.gauss(sigma, infiniteImg, img);
		ImageJFunctions.show(img);
	}
	
	public static void main(String[] args) throws IncompatibleTypeException{
		new ImageJ();
		new Ex6a2();
	}
}
