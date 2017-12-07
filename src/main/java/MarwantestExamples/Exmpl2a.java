package MarwantestExamples;



import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;

public class Exmpl2a {

	public Exmpl2a() throws ImgIOException {

		Img<FloatType> img = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
		Img<FloatType> duplicate = copyImage(img);
		
		ImageJFunctions.show(duplicate);
	}
	
	private <T extends Type<T>> Img<T> copyImage(Img<T> input) {
		Img<T> output = input.factory().create(input, input.firstElement());
		Cursor<T> cursorInput = input.cursor();
		Cursor<T> cursorOutput = output.cursor();
		
		while (cursorInput.hasNext()) {
			cursorInput.fwd();
			cursorOutput.fwd();
			cursorOutput.get().set(cursorInput.get());
		}
		return output;

		
	}

	public static void main(String[] args) throws ImgIOException {
	
		new ImageJ();
		new Exmpl2a();
	}
}
