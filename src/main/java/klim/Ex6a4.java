package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex6a4 {
	public Ex6a4() throws IncompatibleTypeException{
		Img <FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		double[] sigma = new double[img.numDimensions() - 1];
		
		for (int d = 0; d < sigma.length; ++d)
			sigma[d] = 8;
		
		for (int dim = 0; dim < img.numDimensions(); ++dim){
			for (long pos = 0; pos < img.dimension(dim); ++pos){
				if (pos/10 % 2 == 1){
					RandomAccessibleInterval<FloatType> view = Views.hyperSlice(img, dim, pos);
					Gauss3.gauss(sigma, Views.extendMirrorSingle(view), view);
				}
			}
		}
		ImageJFunctions.show(img);
	}
	
	public static void main(String [] args) throws IncompatibleTypeException{
		new ImageJ();
		new Ex6a4();
	}
}
