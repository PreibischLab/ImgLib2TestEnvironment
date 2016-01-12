package stephan;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MinMax
{
	public static < T extends Comparable< T > & Type< T > > RandomAccess< T > max( final RandomAccessibleInterval< T > img )
	{
		final Cursor< T > c = Views.iterable( img ).cursor();
		final RandomAccess< T > max = img.randomAccess();

		// set it to the first element of the img
		c.fwd();

		// take this as first intensity maxima
		max.setPosition( c );

		while( c.hasNext() )
		{
			c.fwd();

			if ( c.get().compareTo( max.get() ) > 0 )
				max.setPosition( c );
		}
		
		return max;
	}

	public static void main( String[] args )
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );
		ImageJFunctions.show( img );
		
		RandomAccess< FloatType > max = max( img );
		
		System.out.println( max.get() );
		System.out.println( Util.printCoordinates( max ) );
		
	}
}
