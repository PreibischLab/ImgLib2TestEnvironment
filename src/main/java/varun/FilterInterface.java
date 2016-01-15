package varun;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class FilterInterface {
	
	
	
	public interface shape {
		
		IterableInterval<FloatType> chooseShape(RandomAccessibleInterval<FloatType> img,Localizable cursor, FinalInterval neighbours);
		
	}
	
	
	
	public static class Rectangle implements shape{
		
		public IterableInterval<FloatType> chooseShape(RandomAccessibleInterval<FloatType> img ,Localizable cursor, FinalInterval neighbours){
			
			final RandomAccessible<FloatType> infinite= Views.extendZero(img);
			
			long min[]=new long[img.numDimensions()];
			
			long max[]=new long[img.numDimensions()];
			
			for (int d=0; d<img.numDimensions(); d++){
				
				min[d]=cursor.getLongPosition(d)-neighbours.dimension(d)/2;
				max[d]=cursor.getLongPosition(d)+neighbours.dimension(d)/2;
				
				
				
			}
			
			
			
			
			FinalInterval interval= new FinalInterval(min,max);
			
			
		final	IterableInterval<FloatType> imgav = Views.interval(infinite, interval) ;
		
		return imgav;
			
		}
		
	}
	
	
	

	public static class Sphere implements shape{
		
		public IterableInterval<FloatType> chooseShape(RandomAccessibleInterval<FloatType> img ,Localizable cursor, FinalInterval neighbours){
			
			final RandomAccessible<FloatType> infinite= Views.extendZero(img);
			
			long radius=0;
			
			long min[]=new long[img.numDimensions()];
			
			long max[]=new long[img.numDimensions()];
			
			for (int d=0; d<img.numDimensions(); d++){
				
				min[d]=cursor.getLongPosition(d)-neighbours.dimension(d)/2;
				max[d]=cursor.getLongPosition(d)+neighbours.dimension(d)/2;
				
				
				
			}
			
			
			for (int d=0; d<img.numDimensions(); ++d){
				radius+=Math.abs(max[d]-min[d])/2;
				
			}
			
			
			
			
		final	HyperSphere<FloatType> imgav=new HyperSphere<FloatType>(infinite, cursor, radius); 
			

		
		return imgav;
			
		}
		
	}
	
	
	public static void meanfilter(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout, final shape rect){
		
		final int n=img.numDimensions();
		final FloatType ini=Views.iterable(img).firstElement().createVariable();
		
		long max[]=new long [n];
		
		long min[]=new long[n];
		
		
		for (int d=0; d<n; ++d){
			
			min[d]=img.min(d);
			
			max[d]=img.max(d);
			
		}
		
		
		
		long maxNeighbours[]= new long[n];
		
		long minNeighbours[]=new long[n];
		
		for(int d=0; d<n; ++d){
			
			
			maxNeighbours[d]=3;
			
			minNeighbours[d]=-3;
		}
		
		
		FinalInterval neighbours= new FinalInterval(minNeighbours,maxNeighbours);
		
		final Cursor<FloatType> bound = Views.iterable(img).cursor();
		
		final RandomAccess<FloatType> outbound= imgout.randomAccess();
		
		
		
		while(bound.hasNext()){
			
			bound.fwd();
		
			final IterableInterval<FloatType> averaged;
			
			averaged=rect.chooseShape(img,bound,neighbours);
			
			
			
			for (FloatType pxval:averaged){
				
				ini.add(pxval);
				
			}
				ini.mul(1.0/averaged.size());

				outbound.setPosition(bound);
				outbound.get().set(ini);
				
			
		}
		
		
		
		
		
		
		
	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		RandomAccessibleInterval<FloatType> imgout=img.copy();

		meanfilter(img,imgout, new Rectangle());
		ImageJFunctions.show(img);
		ImageJFunctions.show(imgout).setTitle("Rectangular_filter");

		meanfilter(img,imgout, new Sphere());

		ImageJFunctions.show(imgout).setTitle("Spherical_filter");

	}

}
