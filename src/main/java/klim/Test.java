package klim;

import java.io.File;

import com.sun.tools.javac.util.Pair;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Test {

	public static Pair<FloatType, FloatType> computeMinMax(RandomAccessibleInterval<FloatType> inputImg){
		final Cursor<FloatType> cursorInput = Views.iterable(inputImg).localizingCursor();
		FloatType type = cursorInput.next();
		FloatType minVal = type.copy();
		FloatType maxVal = type.copy();
		
		while(cursorInput.hasNext()){
			type = cursorInput.next(); 
			if (minVal.get() > type.get())
				minVal.set(type);
			if (maxVal.get() < type.get())
				maxVal.set(type);
				
		}
		
		Pair<FloatType, FloatType> result = new Pair <FloatType, FloatType>(minVal, maxVal);
		
		return result;
		
	}
	
	
	public static void main(String[] args) {
		Img<FloatType> img1 = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		FloatType minVal = new FloatType();
		FloatType maxVal = new FloatType();
		Pair<FloatType, FloatType> result = new Pair <FloatType, FloatType>(minVal, maxVal);
		
		result = computeMinMax(img1);
		System.out.print(result.fst.get() + " " + result.snd.get() + "\n");
		
		ImageJFunctions.show(img1);
	}
}
