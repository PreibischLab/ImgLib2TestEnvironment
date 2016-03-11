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

	/**
	 * Euclidean-Distance Transform based on KD-Tree
	 * 
	 * @param in - binary image
	 * @param out - distances to the closest 1 from any 0
	 */
	public static < T extends RealType< T > > void distanceTransformKD( final IterableInterval< BitType > in, final RandomAccessibleInterval< T > out )
	{
		// make an empty list
		final RealPointSampleList< BitType > list = new RealPointSampleList< BitType >( in.numDimensions() );

		// cursor on the binary image
		final Cursor< BitType > cMain = in.localizingCursor();

		// for every pixel that is 1, make a new RealPoint at that location
		while ( cMain.hasNext() )
			if ( cMain.next().getInteger() == 1 )
				list.add( new RealPoint( cMain ), cMain.get() );

		// build the KD-Tree from the list of points that == 1
		final KDTree< BitType > tree = new KDTree< BitType >( list );

		// Instantiate a nearest neighbor search on the tree (does not modifiy the tree, just uses it)
		final NearestNeighborSearchOnKDTree< BitType > search = new NearestNeighborSearchOnKDTree< BitType >( tree );

		// randomaccess on the output
		final RandomAccess< T > r = out.randomAccess();

		// reset cursor for the input (or make a new one)
		cMain.reset();

		// for every pixel of the binary image
		while ( cMain.hasNext() )
		{
			cMain.fwd();

			// set the randomaccess to the same location
			r.setPosition( cMain );

			// if value == 0, look for the nearest 1-valued pixel
			if ( cMain.get().getInteger() == 0 )
			{
				// search the nearest 1 to the location of the cursor (the current 0)
				search.search( cMain );

				// get the distance (the previous call could return that, this for generality that it is two calls)
				r.get().setReal( search.getDistance() );
			}
			else
			{
				// if value == 1, no need to search because we know the distance is 0
				// but we could search for fun, should return 0
				r.get().setZero();
			}
		}
	}

	public static void main( String[] args ) throws IncompatibleTypeException
	{
		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/test.jpg" ) );

		// take the imgfactory of the img which is of FloatType and transform it into a BitType factory
		final Img< BitType > threshold = img.factory().imgFactory( new BitType() ).create( img, new BitType() );

		Thresholding.threshold( img, threshold, new FloatType( 100 ) );
		ImageJFunctions.show( threshold ).setTitle( "threshold" );

		long t;

		t = System.currentTimeMillis();
		distanceTransformKD( threshold, img );
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
