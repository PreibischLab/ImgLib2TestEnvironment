package marwan;

import java.util.Vector;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;

public class SplitImage {

	public <T extends NumericType<T> & NativeType<T>> SplitImage() throws ImgIOException {
		Helper.log = false;
		Img<FloatType> image = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		ImageJFunctions.show(image);
		
		// We will split the image in 4 columns 2 rows
		final int columns = 4;
		final int rows = 2;
		Vector<RandomAccessibleInterval<FloatType>> views = Helper.splitImage(image, columns, rows);
		for(RandomAccessibleInterval<FloatType> view:views)
			ImageJFunctions.show(view);
	}

	
	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new SplitImage();
	}
}
