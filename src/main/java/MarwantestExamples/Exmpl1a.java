package MarwantestExamples;

import java.io.File;

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;

public class Exmpl1a {

	public <T extends NumericType<T> & NativeType<T>>Exmpl1a() {

		File file  = new File("src/main/resources/DrosophilaWing.tif");
		
		final ImagePlus imp = new Opener().openImage(file.getAbsolutePath());
		
		imp.show();
		
		final Img<T> image = ImagePlusAdapter.wrap(imp);
		
		ImageJFunctions.show(image);
		}
	public static void main(String[] args) {
//		new ImageJ();
		new Exmpl1a();
	}
}
