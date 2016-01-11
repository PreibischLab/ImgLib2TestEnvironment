package ella;

import java.io.File;

import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.Cursor;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;

public class Copy {
	
	
	

	public static void copyImg(final Img< FloatType > img1, final Img< FloatType > img2 ) {
        
		
		final Cursor< FloatType > cursorImg1 = img1.localizingCursor();
        final Cursor< FloatType > cursorImg2 = img2.localizingCursor();
        
        while ( cursorImg1.hasNext())
        {
            cursorImg1.fwd();
            cursorImg2.fwd();
 
            cursorImg2.get().set( cursorImg1.get() );
        }
		
	}
	
	public static void main(String[] args) {
		Img< FloatType > img1 = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		Img<FloatType> img2 = img1.factory().create(img1, img1.firstElement());
		copyImg(img1, img2);
		ImageJFunctions.show(img2);
	}
}
