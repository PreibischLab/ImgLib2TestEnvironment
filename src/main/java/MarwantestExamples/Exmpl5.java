package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.outofbounds.OutOfBoundsFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.ExtendedRandomAccessibleInterval;
import net.imglib2.view.Views;

public class Exmpl5 {

	public Exmpl5() throws ImgIOException {
		Img<FloatType> image = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		
		RandomAccessible<FloatType> infinite1 = Views.extendValue(image, new FloatType(0));
		
		RandomAccessible<FloatType> infinite2 = Views.extendValue(image, new FloatType(128));
		
		RandomAccessible<FloatType> infinite3 = Views.extendRandom(image, 0, 255);
		
		RandomAccessible<FloatType> infinite4 = Views.extendMirrorSingle(image);
		
		RandomAccessible<FloatType> infinite5 = Views.extendMirrorDouble(image);
		
		RandomAccessible<FloatType> infinite6 = Views.extendPeriodic(image);
		
		RandomAccessible<FloatType> infinite7 = new ExtendedRandomAccessibleInterval<FloatType,Img<FloatType>>(image, 
				new OutOfBoundsConstantValueFactory<FloatType,Img<FloatType>>(new FloatType(255)));

		long[] min = new long[image.numDimensions()];
		long[] max = new long[image.numDimensions()];
				
		for (int d = 0; d< image.numDimensions(); ++d) {
			min[d] = -image.dimension(d)-90;
			max[d] = image.dimension(d) * 2 -1 + 90;
			
			FinalInterval interval = new FinalInterval(min,max);

			ImageJFunctions.show(Views.interval(infinite1, interval));
			ImageJFunctions.show(Views.interval(infinite2, interval));
			ImageJFunctions.show(Views.interval(infinite3, interval));
			ImageJFunctions.show(Views.interval(infinite4, interval));
			ImageJFunctions.show(Views.interval(infinite5, interval));
			ImageJFunctions.show(Views.interval(infinite6, interval));
			ImageJFunctions.show(Views.interval(infinite7, interval));
		}
		
	}
	
	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new Exmpl5();
	}
}
