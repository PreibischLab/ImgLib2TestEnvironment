package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.process.ImageProcessor;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.IterableRealInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Sampler;
import net.imglib2.algorithm.neighborhood.*;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetZero;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import net.imglib2.algorithm.region.localneighborhood.*;
import net.imglib2.display.projector.sampler.IntervalSampler;


public class Meanfilter {
	public static Img< FloatType > openAs32Bit( final File file )
	{
		return openAs32Bit( file, new ArrayImgFactory< FloatType >() );
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

	
	public static  void mean ( RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout){
	
	//	final RandomAccessibleInterval <FloatType> imgout;
		
		
		Interval interval = Intervals.expand( img, -100 );
		
		 img=Views.interval(img, interval);
		 
		
			
		 imgout=Views.interval(imgout, interval);
		 
		 final Cursor < FloatType > bound = Views.iterable(img).cursor();
		 
		 final Cursor < FloatType > outbound = Views.iterable(imgout).cursor();
		
		// FloatType mean=((IterableRealInterval<FloatType>) img).firstElement();
		 FloatType mean= new FloatType();
		 
		 mean.setZero();
		 
		 
		
	
		 
		// value.setReal(ini);
		 
		 while(bound.hasNext()){
			 
		  
		  bound.fwd();
		  outbound.fwd();
		  
		  mean.add(bound.get());
         
		 }
		 
		 outbound.get().set(mean);
		
		
		
	
	}
	
	
	public static void main (String[]args){
		
		final Img< FloatType > img = openAs32Bit( new File( "src/main/resources/bridge.png" ) );
	final Img < FloatType > imgout=img.copy();
		
		mean(img,imgout);
		ImageJFunctions.show(img);
		ImageJFunctions.show(imgout);
	}

}
