package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.RandomAccess;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class Ex2b {

	public <T extends Type <T>> Img <T> copyImageWrong (final Img< T > input,
			final ImgFactory< T > imgFactory ){
		Img <T> output = imgFactory.create(input, input.firstElement());
		
		Cursor <T> cursorInput = input.cursor();
		Cursor <T> cursorOutput = output.cursor();
		
		while (cursorInput.hasNext()){
			cursorInput.fwd();
			cursorOutput.fwd();
			
			cursorOutput.get().set(cursorInput.get());
		}
		
		return output;
	}
	
	public <T extends Type <T>> Img <T> copyImageRight(final Img< T > input,
			final ImgFactory< T > imgFactory ){
		Img <T> output = imgFactory.create(input, input.firstElement());
		Cursor <T> cursorInput = input.cursor();
		RandomAccess <T> randomAccess = output.randomAccess();
		while (cursorInput.hasNext()){
			cursorInput.fwd();
			randomAccess.setPosition(cursorInput);
			randomAccess.get().set(cursorInput.get());
		}		
		
		return output;
		
	}
	public  Ex2b(){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		Img<FloatType> imgOut1 = copyImageWrong(img, new CellImgFactory< FloatType >( 20 ) );
		Img<FloatType> imgOut2 = copyImageRight(img, new CellImgFactory< FloatType >( 20 ) );
		ImageJFunctions.show(img).setTitle("initial");	
		ImageJFunctions.show(imgOut1).setTitle("wrong");
		ImageJFunctions.show(imgOut2).setTitle("right");
	} 
	
	public static void main (String[] args){
		new ImageJ();
		new Ex2b();
		
		System.out.println("We\'re are done!");
		
	}
	
	
}
