package david;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImg;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class DavidTest {

	public static <T extends Type<T>> Img<T> copyImage(Img<T> source){
		
		Img<T> dest = source.factory().create(source, source.firstElement());
		Cursor<T> destC = dest.cursor();
		Cursor<T> srcC = source.cursor();
		
		
		while (destC.hasNext()){
			destC.fwd();
			srcC.fwd();
			destC.get().set(srcC.get());
		}
		
		return dest;
		
	}
	
	public static <T extends Type<T>> void copyImage(RandomAccessibleInterval<T> source, IterableInterval<T> dest) {
		Cursor<T> destC = dest.localizingCursor();
		RandomAccess<T> srcRA = source.randomAccess();
		
		
		while (destC.hasNext()){
			destC.fwd();
			srcRA.setPosition(destC);;
			destC.get().set(srcRA.get());
		}
	}
	
	public static void main(String[] args) {
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		Img<FloatType> img2 = new CellImgFactory<FloatType>().create(img, new FloatType());
		ImageJFunctions.show(img2);
		copyImage(img, img2);
		ImageJFunctions.show(img2);
		
	}
	
	

}
