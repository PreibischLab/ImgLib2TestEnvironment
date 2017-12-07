package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Exmpl6a1 {

	public Exmpl6a1() throws ImgIOException, IncompatibleTypeException {
		Img<FloatType> image = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		ImageJFunctions.show(image);
		double[] sigma = new double [ image.numDimensions()];
		
		for (int d = 0; d< image.numDimensions();++d) {
			sigma[d]=8;
			
		}
		RandomAccessible<FloatType> infiniteImg = Views.extendValue(image, new FloatType());
		
		Gauss3.gauss(sigma, infiniteImg, image);
//		ImageJFunctions.show(Gauss.toFloat(sigma, image));

		ImageJFunctions.show( image);
		
	}
	public static void main(String[] args) throws ImgIOException, IncompatibleTypeException {
		new ImageJ();
		new Exmpl6a1();
		
	}
}
