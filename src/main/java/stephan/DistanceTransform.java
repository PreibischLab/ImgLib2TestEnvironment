package stephan;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

/**
 * Based on a binary image (from thresholding), compute for every bright pixel the closest distance to a 0-pixel 
 * 
 *
 */
public class DistanceTransform
{
	public static double dist( final Localizable l1, final Localizable l2 )
	{
		double dist = 0;

		for ( int d = 0; d < l1.numDimensions(); ++d )
			dist += Math.pow( l2.getDoublePosition( d ) - l1.getDoublePosition( d ), 2 );

		return Math.sqrt( dist );
	}

	public static < T extends RealType< T > > void distanceTransform( final IterableInterval< BitType > in, final RandomAccessibleInterval< T > out )
	{
		final Cursor< BitType > cMain = in.localizingCursor();
		final RandomAccess< T > r = out.randomAccess();

		while ( cMain.hasNext() )
		{
			cMain.fwd();
			r.setPosition( cMain );

			// if this pixel is "0" iterate over all to see which pixel with value "1" is closest.
			if ( cMain.get().getInteger() == 0 )
			{
				final Cursor< BitType > cLocal = in.cursor();
				double minDist = Double.MAX_VALUE;

				while( cLocal.hasNext() )
					if ( cLocal.next().getInteger() == 1 )
						minDist = Math.min( minDist, dist( cLocal, cMain ) );
	
				r.get().setReal( minDist );
			}
			else
			{
				r.get().setReal( 0 );
			}
		}
	}

	public static void main( String[] args ) throws IncompatibleTypeException
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/dt.png" ) );

		// take the imgfactory of the img which is of FloatType and transform it into a BitType factory
		final Img< BitType > threshold = img.factory().imgFactory( new BitType() ).create( img, new BitType() );

		Thresholding.threshold( img, threshold, new FloatType( 200 ) );
		ImageJFunctions.show( threshold );

		distanceTransform( threshold, img );
		ImageJFunctions.show( img );
	}
}
