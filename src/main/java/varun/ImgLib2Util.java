package varun;

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class ImgLib2Util
{
	public static Img< FloatType > openAs32Bit( final File file )
	{
		return openAs32Bit( file, new ArrayImgFactory< FloatType >() );
	}

public static Img<FloatType> Copy32Bit(final Img< FloatType > img1)
{

final Img< FloatType > img2;


			img2 = img1.factory().create( img1, new FloatType ());


final Cursor< FloatType > c1 = img1.localizingCursor();

final Cursor< FloatType > c2 = img2.localizingCursor();

while ( c1.hasNext() )
			{

c1.fwd();
c2.fwd();


 c2.get().set(c1.get());




}

return img2;

}


	public static Img< FloatType > openAs32Bit( final File file, final ImgFactory< FloatType > factory )
	{
		if ( !file.exists() )
			throw new RuntimeException( "File '" + file.getAbsolutePath() + "' does not exisit." );

		final ImagePlus imp = new Opener().openImage( file.getAbsolutePath() );

		if ( imp == null )
			throw new RuntimeException( "File '" + file.getAbsolutePath() + "' coult not be opened." );

		final Img< FloatType > img;

		if ( imp.getStack().getSize() == 1 )
		{
			// 2d
			img = factory.create( new int[]{ imp.getWidth(), imp.getHeight() }, new FloatType() );
			final ImageProcessor ip = imp.getProcessor();

			final Cursor< FloatType > c = img.localizingCursor();
			
			while ( c.hasNext() )
			{
				c.fwd();

				final int x = c.getIntPosition( 0 );
				final int y = c.getIntPosition( 1 );

				c.get().set( ip.getf( x, y ) );
			}

		}
		else
		{
			// >2d
			img = factory.create( new int[]{ imp.getWidth(), imp.getHeight(), imp.getStack().getSize() }, new FloatType() );

			final Cursor< FloatType > c = img.localizingCursor();

			// for efficiency reasons
			final ArrayList< ImageProcessor > ips = new ArrayList< ImageProcessor >();

			for ( int z = 0; z < imp.getStack().getSize(); ++z )
				ips.add( imp.getStack().getProcessor( z + 1 ) );

			while ( c.hasNext() )
			{
				c.fwd();

				final int x = c.getIntPosition( 0 );
				final int y = c.getIntPosition( 1 );
				final int z = c.getIntPosition( 2 );

				c.get().set( ips.get( z ).getf( x, y ) );
			}
		}

		return img;
	}

	public static void main( String[] args )
	{
		new ImageJ();
		
		//final Img< FloatType > img = openAs32Bit( new File( "src/main/resources/mri-stack.tif" ) );
		final Img< FloatType > img = openAs32Bit( new File( "src/main/resources/bridge.png" ) );
		
		ImageJFunctions.show( img );


final Img< FloatType > img2=Copy32Bit(img);

ImageJFunctions.show(img2);

	}
}
