package varun;

import java.io.File;



import net.imglib2.Cursor;
import net.imglib2.FinalInterval;

import net.imglib2.IterableInterval;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;

import net.imglib2.img.Img;

import net.imglib2.img.display.imagej.ImageJFunctions;

import net.imglib2.type.numeric.real.FloatType;


import net.imglib2.view.Views;

public class Meanfilter {
	

	public static void meanImage(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval< FloatType > imgout) {

		
		
		final RandomAccessible<FloatType> infinite = Views.extendZero(img);

		
		final FloatType ini=Views.iterable(img).firstElement().createVariable();
		
		
		
		
		final int n = img.numDimensions();
		long min[] = new long[n];
		long max[] = new long[n];

		long nearmin[]= new long [n];
		long nearmax[]=new long[n];
		
		for(int d=0; d<n; ++d){
			 
			nearmin[d]=-3;
			
			nearmax[d]=3;
			
		}
		
		final FinalInterval neighbours= new FinalInterval(nearmin,nearmax); // Determines the nearest neighbors.
		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();
		final RandomAccess< FloatType > outbound=imgout.randomAccess();
	
		
		
		
		while(bound.hasNext()){
		
			bound.fwd();
			
		for (int d = 0; d < n; ++d) {
			min[d] = bound.getLongPosition(d)-neighbours.dimension(d)/2;
			max[d] = bound.getLongPosition(d)+neighbours.dimension(d)/2;

		}

		

		FinalInterval interval = new FinalInterval(min, max);
		
	final IterableInterval<FloatType>	imgav = Views.interval(infinite, interval);
	
	
	//final HyperSphere< FloatType > imgav= new HyperSphere<FloatType>(infinite, bound, 3);

	for (FloatType pxval:imgav){
		
		ini.add(pxval);
		
	}
		ini.mul(1.0/imgav.size());

		outbound.setPosition(bound);
		outbound.get().set(ini);
		
		
		}
		
	
		
		
	}

	

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		RandomAccessibleInterval<FloatType> imgout=img.copy();

		meanImage(img,imgout);
		ImageJFunctions.show(img);
		ImageJFunctions.show(imgout);

		

		//computeNeighbourhood(img, imgout);

	}

}
