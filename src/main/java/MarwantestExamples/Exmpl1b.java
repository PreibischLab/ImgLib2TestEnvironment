package MarwantestExamples;

import java.io.File;

import ij.ImageJ;
import io.scif.config.SCIFIOConfig;
import io.scif.config.SCIFIOConfig.ImgMode;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class Exmpl1b {

	public <T extends RealType<T> & NativeType<T>> Exmpl1b() throws ImgIOException {
		File file = new File("src/main/resources/DrosophilaWing.tif");
		String path = file.getAbsolutePath();
		ImgOpener imgOpener = new ImgOpener();
		Img<T> image = (Img<T>) imgOpener.openImg(path);
		ImageJFunctions.show(image);
		SCIFIOConfig config = new SCIFIOConfig();
		config.imgOpenerSetImgModes(ImgMode.CELL);
		
		Img<T> imageCll = (Img<T>) imgOpener.openImg(path, config);
		
		ImageJFunctions.show(imageCll);
	}
	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new Exmpl1b();
	}
}
