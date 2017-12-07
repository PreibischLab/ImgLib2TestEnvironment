package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;

public class Exmpl2b {
 
	public Exmpl2b() throws ImgIOException {

		Img<FloatType> img = (Img<FloatType>) new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new ArrayImgFactory<FloatType>());
	
Img<FloatType> duplicate = copyImageCorrect(img, new CellImgFactory<FloatType>(20));
ImageJFunctions.show(img);
		ImageJFunctions.show(duplicate);
	
	}
	private <T extends Type<T>> Img<T> copyImageCorrect(Img<T> input, ImgFactory<T> imgFactory) {
		
		Img<T> output = input.factory().create(input, input.firstElement());
		Cursor<T> cursorInput = input.localizingCursor();
		RandomAccess<T> randomAccess = output.randomAccess();
		
		while (cursorInput.hasNext()) {
			cursorInput.fwd();
			randomAccess.setPosition(cursorInput);;
			randomAccess.get().set(cursorInput.get());
		}
		return output;
		
	}
	public static void main(String[] args) throws ImgIOException {
	new ImageJ();
	new Exmpl2b();
	}
}
