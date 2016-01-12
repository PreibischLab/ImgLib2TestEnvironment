package stephan;

import java.util.ArrayList;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class Sphere
{
	public static < T extends Type< T > & Comparable< T > > void drawSphere(
			final RandomAccessible< T > img,
			final T value,
			final Localizable center,
			final long radius )
	{
		for ( final T t : new HyperSphere< T >( img, center, radius ) )
			if ( value.compareTo( t ) > 0 )
				t.set( value );
	}


	public static void main( String[] args )
	{
		new ImageJ();
		final Img< FloatType > img = ArrayImgs.floats( 512, 512 );

		final Cursor< FloatType > c = new HyperSphere< FloatType >( img, new Point( 256, 256 ), 200 ).localizingCursor();

		final Random rnd = new Random( 5343 );

		while ( c.hasNext() )
		{
			c.fwd();

			if ( rnd.nextDouble() > 0.98 )
				drawSphere(
						img,
						new FloatType( rnd.nextFloat() ),
						c,
						Math.round( rnd.nextDouble() * 10 )  );
		}

		ImageJFunctions.show( img );
	}
}
