package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class Exmpl6a2 {

	public Exmpl6a2() throws ImgIOException, IncompatibleTypeException {
		Img<FloatType> image = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		
		double[] sigma = new double [ image.numDimensions()];
		
		for (int d = 0; d< image.numDimensions();++d) {
			System.out.println(d+"|");
			sigma[d]=0;
		}
		RandomAccessible<FloatType> infiniteImg = Views.extendMirrorSingle(image);
		long[] min = new long[image.numDimensions()];
		long[] max = new long[image.numDimensions()];
				
//		for (int d = 0; d< image.numDimensions(); ++d) {
//			min[d] = -image.dimension(d)-90;
//			max[d] = image.dimension(d) * 2 -1 + 90;
//			
//			FinalInterval interval = new FinalInterval(min,max); 
			

		FinalInterval interval = Intervals.createMinMax(300,30,500,250);
		ImageJFunctions.show(Views.interval(infiniteImg, interval));
		RandomAccessibleInterval<FloatType> region = Views.interval(image, interval);
		
		Gauss3.gauss(sigma, infiniteImg, region);
//		ImageJFunctions.show(Gauss.toFloat(sigma, image));

		ImageJFunctions.show( image);
		
	}
	public static void main(String[] args) throws ImgIOException, IncompatibleTypeException {
		new ImageJ();
		new Exmpl6a2();
		
	}
}
