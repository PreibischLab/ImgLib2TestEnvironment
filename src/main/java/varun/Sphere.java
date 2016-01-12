package varun;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class Sphere
{

	public static < T extends RealType < T >> void drawMulti( 
			final RandomAccessibleInterval < T > img, Point newcenter, long newradius)
	{
		
		
		HyperSphere< T > bigSphere =
                new HyperSphere< T >( img, newcenter, newradius );
		 
		  
	    

		Cursor< T > bigc = bigSphere.cursor();

		// long smvar = 0;
		
		long newradiustwo = newradius/10;
		while ( bigc.hasNext() )
		{
			
			
			bigc.fwd();
			
			
			//smallc.get().setReal( smvar );
			HyperSphere< T > smallSphere =
	                new HyperSphere< T >( img, bigc, newradiustwo );
			
			//++smvar;
			
			
			
			

}
		
	
		
		
	}
	
	public static < T extends RealType< T >> void drawSphere(
			final RandomAccessibleInterval< T > img )
	{

		// the number of dimensions
		long[] center = new long[ img.numDimensions() ];
		for (int d = 0; d < img.numDimensions(); ++d)
		{
			center[ d ] = ( img.max( d ) - img.min( d ) ) / 2 + img.min( d );
		}

		// define the center and radius
		Point centerdef = new Point( center );
		long radius = 10;
		HyperSphere< T > testsphere = new HyperSphere< T >( img, centerdef,
				radius );

		Cursor< T > c = testsphere.cursor();

		long var = 0;

		while ( c.hasNext() )
		{
			c.fwd();
			c.get().setReal( var );
			++var;

		}

	}

	public static void main( String[] args )
	{
		new ImageJ();

		final Img< FloatType > img = ArrayImgs.floats( 100, 100 );
		drawSphere( img );
		ImageJFunctions.show( img );
		
		Point newcenter = new Point( 50 );
		drawMulti(img, newcenter,2);
		
	}

}
