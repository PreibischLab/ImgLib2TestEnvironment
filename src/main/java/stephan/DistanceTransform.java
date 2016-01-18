package stephan;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
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
	public interface Distance
	{
		double dist( final Localizable l1, final Localizable l2 );
	}

	/**
	 * l1-norm
	 */
	public static class ManhattanDistance implements Distance
	{
		@Override
		public double dist( Localizable l1, Localizable l2 )
		{
			double dist = 0;

			for ( int d = 0; d < l1.numDimensions(); ++d )
				dist += Math.abs( l2.getDoublePosition( d ) - l1.getDoublePosition( d ) );

			return dist;
		}
	}

	/**
	 * l2-norm
	 */
	public static class EuclideanDistance implements Distance
	{
		@Override
		public double dist( Localizable l1, Localizable l2 )
		{
			double dist = 0;

			for ( int d = 0; d < l1.numDimensions(); ++d )
				dist += Math.pow( l2.getDoublePosition( d ) - l1.getDoublePosition( d ), 2 );

			return Math.sqrt( dist );
		}
	}

	/**
	 * l_infinite?-norm / max-distance
	 */
	public static class ChessboardDistance implements Distance
	{
		@Override
		public double dist( Localizable l1, Localizable l2 )
		{
			double dist = Math.abs( l2.getDoublePosition( 0 ) - l1.getDoublePosition( 0 ) );

			for ( int d = 1; d < l1.numDimensions(); ++d )
				dist = Math.max( dist, Math.abs( l2.getDoublePosition( d ) - l1.getDoublePosition( d ) ) );

			return dist;
		}
	}

	public static < T extends RealType< T > > void distanceTransform( final IterableInterval< BitType > in, final RandomAccessibleInterval< T > out, final Distance distance )
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
						minDist = Math.min( minDist, distance.dist( cLocal, cMain ) );
	
				r.get().setReal( minDist );
			}
			else
			{
				r.get().setReal( 0 );
			}
		}
	}

	public static < T extends RealType< T > > void distanceTransformKD( final IterableInterval< BitType > in, final RandomAccessibleInterval< T > out, final Distance distance )
	{
		final RealPointSampleList< BitType > list = new RealPointSampleList< BitType >( in.numDimensions() );
		
		final Cursor< BitType > cMain = in.localizingCursor();

		while ( cMain.hasNext() )
			if ( cMain.next().getInteger() == 1 )
				list.add( new RealPoint( cMain ), cMain.get() );

		final KDTree< BitType > tree = new KDTree< BitType >( list );
		final NearestNeighborSearchOnKDTree< BitType > search = new NearestNeighborSearchOnKDTree< BitType >( tree );

		final RandomAccess< T > r = out.randomAccess();
		cMain.reset();

		while ( cMain.hasNext() )
		{
			cMain.fwd();
			r.setPosition( cMain );

			if ( cMain.get().getInteger() == 0 )
			{
				search.search( cMain );
				r.get().setReal( search.getDistance() );
			}
			else
			{
				r.get().setZero();
			}
		}
	}

	public static void main( String[] args ) throws IncompatibleTypeException
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );

		// take the imgfactory of the img which is of FloatType and transform it into a BitType factory
		final Img< BitType > threshold = img.factory().imgFactory( new BitType() ).create( img, new BitType() );

		Thresholding.threshold( img, threshold, new FloatType( 200 ) );
		ImageJFunctions.show( threshold ).setTitle( "threshold" );

		long t;

		t = System.currentTimeMillis();
		distanceTransformKD( threshold, img, new EuclideanDistance() );
		System.out.println( "O(n logn): " + (System.currentTimeMillis() - t) + " ms." );
		ImageJFunctions.show( img ).setTitle( "kd-euclidean distance" );

		t = System.currentTimeMillis();
		distanceTransform( threshold, img, new EuclideanDistance() );
		System.out.println( "O(n^2): " + (System.currentTimeMillis() - t) + " ms." );
		ImageJFunctions.show( img ).setTitle( "euclidean distance" );


		/*
		distanceTransform( threshold, img, new ManhattanDistance() );
		ImageJFunctions.show( img ).setTitle( "manhattan distance" );

		distanceTransform( threshold, img, new ChessboardDistance() );
		ImageJFunctions.show( img ).setTitle( "chessboard distance" );*/
	}
}
