package ella;

import java.io.File;

import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;

public class Copy {
	
<<<<<<< HEAD
=======
	
	
	public static <T extends Type<T>> Img< T > copyImg(final Img<T> img1 )
	{
		final Img< T > img2 = img1.factory().create( img1, img1.firstElement() );
		copyImg( img1, img2 );
		return img2;
	}
>>>>>>> 888cce051aa3e004b6ff61531267bd522d2e5343

	public static <T extends Type<T>> void copyImg(final IterableInterval<T> img1, final RandomAccessibleInterval<T> img2 )
	{

		if ( img1.iterationOrder().equals( Views.iterable( img2 ).iterationOrder() ))
		{
			final Cursor< T > cursorImg1 = img1.cursor();
	        final Cursor< T > cursorImg2 = Views.iterable( img2 ).cursor();
	        
	        while ( cursorImg1.hasNext())
	        {
	            cursorImg1.fwd();
	            cursorImg2.fwd();
	 
	            cursorImg2.get().set( cursorImg1.get() );
	        }
		}
		else
		{
			final Cursor< T > cursorImg1 = img1.localizingCursor();
			final RandomAccess< T > ra2 = img2.randomAccess();

			//final int[] pos = new int[ img1.numDimensions() ];

	        while ( cursorImg1.hasNext())
	        {
	            cursorImg1.fwd();
	            ra2.setPosition( cursorImg1 );
	            //cursorImg1.localize( pos );
	            //ra2.setPosition( pos );

	            ra2.get().set( cursorImg1.get() );
	        }
		}
	}
	
	public static void main(String[] args) {
		Img<FloatType> img1 = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		Img<FloatType> img2 = new CellImgFactory< FloatType >().create(img1, img1.firstElement());
		
		copyImg(img1, img2);

		ImageJFunctions.show(  Views.zeroMin( Views.interval( img1, new long[]{ 30,  40 }, new long[]{ 100, 100 } ) ) );
		
		copyImg(
				Views.zeroMin( Views.interval( img1, new long[]{ 30, 40 }, new long[]{ 100, 100 } ) ),
				Views.zeroMin( Views.interval( img2, new long[]{ 330, 140 }, new long[]{ 400, 200 } ) ) );
		
		ImageJFunctions.show(img2);

	}
}
