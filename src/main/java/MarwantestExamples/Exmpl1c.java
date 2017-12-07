package MarwantestExamples;

import ij.ImageJ;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class Exmpl1c {

	public Exmpl1c() {

		final ImgFactory<FloatType> imgFactory = new CellImgFactory<FloatType>(5);
		
		final Img<FloatType> img1 = imgFactory.create(new long[] {20,30,40}, new FloatType());
		final Img<FloatType> img2 = imgFactory.create(img1, img1.firstElement());
		ImageJFunctions.show(img1);
		ImageJFunctions.show(img2);
		
	}
	public static void main(String[] args) {
		new ImageJ();
		new Exmpl1c();
	}
}
