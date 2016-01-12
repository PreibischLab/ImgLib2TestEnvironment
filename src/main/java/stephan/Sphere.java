package stephan;

import java.util.ArrayList;

import ij.ImageJ;
import net.imglib2.Point;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class Sphere
{
	public static void main( String[] args )
	{
		new ImageJ();
		final Img< FloatType > img = ArrayImgs.floats( 100, 100 );

		final HyperSphere< FloatType > sphere = new HyperSphere< FloatType >( img, new Point( 50, 50 ), 10 );

		int i = 0;

		for ( final FloatType t : sphere )
			t.set( i++ );

		ImageJFunctions.show( img );
	}
}
