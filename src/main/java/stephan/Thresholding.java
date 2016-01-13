package stephan;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

/**
 * 
 * For every pixel see if it is above or below a threshold
 * 
 * if above, set 1, otherwise 0 > use a binary image for this (BitType)
 *
 */
public class Thresholding
{
	public static < T extends Comparable< T > > void threshold( final RandomAccessibleInterval< T > in, final RandomAccessibleInterval< BitType > out, final T threshold )
	{
		final Cursor< T > cIn = Views.iterable( in ).localizingCursor();
		final RandomAccess< BitType > rOut = out.randomAccess();

		while ( cIn.hasNext() )
		{
			cIn.fwd();
			rOut.setPosition( cIn );

			if ( cIn.get().compareTo( threshold ) > 0 )
				rOut.get().setOne();
			else
				rOut.get().setZero();
		}
	}
	
	public static void main( String[] args ) throws IncompatibleTypeException
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );

		// take the imgfactory of the img which is of FloatType and transform it into a BitType factory
		final Img< BitType > threshold = img.factory().imgFactory( new BitType() ).create( img, new BitType() );

		threshold( img, threshold, new FloatType( 200 ) );
		ImageJFunctions.show( threshold );
	}
}
