package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import ucar.unidata.geoloc.projection.FlatEarth;

public class Exmpl1d {

	public Exmpl1d() throws ImgIOException {

		Img<FloatType> img = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		ImageJFunctions.show(img);
		RandomAccessibleInterval<FloatType> view  = Views.interval(img, new long[] {200,200}, new long[] {500,350});
		ImageJFunctions.show(view);
		ImageJFunctions.show(Views.rotate(view, 0, 1));
		
	}

	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new Exmpl1d();
	}
}
