package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Exmpl2c {
 
	public Exmpl2c() throws ImgIOException
    {
        // open with ImgOpener as a float
        Img<FloatType> img = new ImgOpener().openImg("DrosophilaWing.tif",
            new FloatType());
 
        // copy & display an image
        Img< FloatType > duplicate = img.factory().create( img, img.firstElement() );
        copy( img, duplicate );
        ImageJFunctions.show( duplicate );
 
        // use a View to define an interval as source for copying
        //
        // Views.offsetInterval() does not only define where it is, but also adds a translation
        // so that the minimal coordinate (upper left) of the view maps to (0,0)
        RandomAccessibleInterval< FloatType > viewSource = Views.offsetInterval( img,
            new long[] { 100, 100 }, new long[]{ 250, 150 } );
 
        // and as target
        RandomAccessibleInterval< FloatType > viewTarget = Views.offsetInterval( img,
            new long[] { 500, 200 }, new long[]{ 250, 150 } );
 
        // now we make the target iterable
        // (which is possible because it is a RandomAccessibleInterval)
        IterableInterval< FloatType > iterableTarget = Views.iterable( viewTarget );
 
        // copy it into the original image (overwriting part of img)
        copy( viewSource, iterableTarget );
 
        // show the original image
        ImageJFunctions.show( img );
    }
	  public < T extends Type< T > > void copy( final RandomAccessible< T > source,
		        final IterableInterval< T > target )
		    {
		        // create a cursor that automatically localizes itself on every move
		        Cursor< T > targetCursor = target.localizingCursor();
		        RandomAccess< T > sourceRandomAccess = source.randomAccess();
		 
		        // iterate over the input cursor
		        while ( targetCursor.hasNext())
		        {
		            // move input cursor forward
		            targetCursor.fwd();
		 
		            // set the output cursor to the position of the input cursor
		            sourceRandomAccess.setPosition( targetCursor );
		 
		            // set the value of this pixel of the output image, every Type supports T.set( T type )
		            targetCursor.get().set( sourceRandomAccess.get() );
		        }
		    }
	public static void main(String[] args) throws ImgIOException {
	new ImageJ();
	new Exmpl2c();
	}
}
