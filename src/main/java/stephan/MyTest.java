package stephan;

import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class MyTest
{
	public static void iterateAllPixels()
	{
		Img< FloatType > test = new CellImgFactory< FloatType >( 3 ).create( new int[]{ 10, 10 }, new FloatType() );
		
		final Cursor< FloatType > c = test.cursor();
		
		int i = 0;
		
		while (  c.hasNext() )
		{
			c.fwd();
			c.get().set( i );
			++i;
		}
		
		ImageJFunctions.show( test );

	}
}
