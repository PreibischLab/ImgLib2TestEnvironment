package stephan;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

/**
 * For every pixel of an image, compute the mean of its neighboring pixels
 * 
 * 
 */
public class MeanFilter
{
	/**
	 * Compute the mean of an iterable interval
	 * 
	 * @param toFilter
	 * @param result
	 */
	public static < T extends RealType<  T > > void filter( final IterableInterval< T > toFilter, final T result )
	{
		result.setZero();

		for ( final T pixelValue : toFilter )
			result.add( pixelValue );

		result.mul( 1.0 / toFilter.size() );
	}

	/**
	 * filter the whole randomaccessible with a kernel defined by an interval
	 * 
	 * @param in
	 * @param out
	 * @param kernel
	 */
	public static < T extends RealType<  T > > void filterRect( final RandomAccessibleInterval< T > in, final RandomAccessibleInterval< T > out, final Interval kernel )
	{
		final RandomAccessible< T > infinite = Views.extendMirrorSingle( in );

		final Cursor< T > c = Views.iterable( in ).localizingCursor();
		final RandomAccess<  T > r = out.randomAccess();

		final int n = in.numDimensions();

		final long[] min = new long[ n ];
		final long[] max = new long[ n ];

		final T tmp = Views.iterable( in ).firstElement().createVariable();

		while ( c.hasNext() )
		{
			c.fwd();

			for ( int d = 0; d < n; ++d )
			{
				min[ d ] = c.getLongPosition( d ) - kernel.dimension( d )/2;
				max[ d ] = c.getLongPosition( d ) + kernel.dimension( d )/2;
			}

			final IterableInterval< T > toAvg = Views.interval( infinite, min, max );

			filter( toAvg, tmp );

			r.setPosition( c );
			r.get().set( tmp );
		}
	}

	/**
	 * filter the whole randomaccessible with a kernel defined by a sphere
	 * 
	 * @param in
	 * @param out
	 * @param kernel
	 */
	public static < T extends RealType<  T > > void filterCircle( final RandomAccessibleInterval< T > in, final RandomAccessibleInterval< T > out, final long radius )
	{
		final RandomAccessible< T > infinite = Views.extendMirrorSingle( in );

		final Cursor< T > c = Views.iterable( in ).localizingCursor();
		final RandomAccess<  T > r = out.randomAccess();

		final T tmp = Views.iterable( in ).firstElement().createVariable();

		while ( c.hasNext() )
		{
			c.fwd();

			final HyperSphere< T > sphere = new HyperSphere< T >( infinite, c, radius );

			filter( sphere, tmp );

			r.setPosition( c );
			r.get().set( tmp );
		}
	}

	public static void main( String[] args )
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );
		final Img< FloatType > outRect = img.copy();
		final Img< FloatType > outSphere = img.copy();

		final int n = img.numDimensions();
		final long[] min = new long[ n ];
		final long[] max = new long[ n ];

		for ( int d = 0; d < n; ++d )
		{
			min[ d ] = -3;
			max[ d ] = 3;
		}

		ImageJFunctions.show( img );

		filterRect( img, outRect, new FinalInterval( min, max ) );
		ImageJFunctions.show( outRect );

		filterCircle( img, outSphere, 3 );
		ImageJFunctions.show( outSphere );
	}
}
