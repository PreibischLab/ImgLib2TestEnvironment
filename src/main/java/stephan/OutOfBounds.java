package stephan;

import java.io.File;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class OutOfBounds
{
	public static void main( String[] args )
	{

		// defined in an interval (RandomAccessibleInterval)
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );

		ImageJFunctions.show( img );

		// transform it into an RandomAccessible
		final RandomAccessible< FloatType > infinite = Views.extendRandom( img, 0, 255 );

		final int n = img.numDimensions();
		long[] min = new long[ n ];
		long[] max = new long[ n ];

		for ( int d = 0; d < n; ++d )
		{
			min[ d ] = img.min( d ) - ( img.dimension( d ) - 1 )/2;
			max[ d ] = img.max( d ) + ( img.dimension( d ) - 1 )/2;
		}

		RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
		//interval = Views.interval( infinite, Intervals.expand( img, -50 ) );
		
		ImageJFunctions.show( interval );
	}
}
