package varun;

import java.io.File;

import net.imglib2.FinalInterval;

import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class zoom {
	
	
	public static RandomAccessibleInterval<FloatType> getZoom (RandomAccessibleInterval<FloatType> img) {
		
		
		final RandomAccessible<FloatType> infinite= Views.extendZero(img);
		
		final int n=img.numDimensions();
		
		long min[]=new long[n];
		long max[]= new long[n];
		
		long shortmin[]=new long [n];
		
		long shortmax[]= new long [n];
		
for( int d=0; d<n; ++d) {
			
			min[d]=img.min(d);
			max[d]=img.max(d);
		}
		
		
		
		for( int d=0; d<n; ++d) {
			
			shortmin[d]=-100;
			shortmax[d]=100;
		}
		
		
		FinalInterval fullinterval= new FinalInterval(min,max);
		
		FinalInterval interval= new FinalInterval(shortmin,shortmax);
		
	
	//	final RandomAccessibleInterval<FloatType> imgout= Views.interval(infinite,fullinterval);
		
	final RandomAccessibleInterval<FloatType>	imgout=Views.interval(infinite, interval);
	
	
	final RandomAccessible<FloatType> tmpimg=Views.extendZero(imgout);
	
	final RandomAccessibleInterval<FloatType> imgzoom= Views.interval(tmpimg, fullinterval);
	
		System.out.println(imgout.dimension(0));
		
		System.out.println(imgzoom.dimension(0));
	return imgzoom;
		
		
		
		
	}
	
	
	
	
	public static void main(String[] args) {
		
		final Img < FloatType > image = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		
		RandomAccessibleInterval<FloatType> imgout=image.copy();
		
		imgout=getZoom(image);
		
		ImageJFunctions.show(imgout);
	}
	

}
