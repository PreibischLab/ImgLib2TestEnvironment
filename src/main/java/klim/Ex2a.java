package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class Ex2a {

	public <T extends Type <T>> Img <T> copyImage (final Img <T> input){
		Img <T> output = input.factory().create(input, input.firstElement());
		
		Cursor <T> cursorInput = input.cursor();
		Cursor <T> cursorOutput = output.cursor();
		
		while (cursorInput.hasNext()){
			cursorInput.fwd();
			cursorOutput.fwd();
			
			cursorOutput.get().set(cursorInput.get());
		}
		
		return output;
	}
	
	public  Ex2a(){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		Img<FloatType> duplicate = copyImage(img);
		ImageJFunctions.show(duplicate);
	} 
	
	public static void main (String[] args){
		new ImageJ();
		new Ex2a();
		
		System.out.println("We\'re are done!");
		
	}
	
	
}
