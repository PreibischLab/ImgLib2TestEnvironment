package klim;

import java.io.File;

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;

import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;

// simple example to show an image 

public class Ex1a {
	
	public < T extends NumericType< T > & NativeType< T >> Ex1a(){
		File file = new File("src/main/resources/test.jpg");
		final ImagePlus imp = new Opener().openImage(file.getAbsolutePath());
		imp.show();
		final Img < T > image = ImagePlusAdapter.wrap(imp);
		ImageJFunctions.show( image );
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Ex1a();
	}
	
}
