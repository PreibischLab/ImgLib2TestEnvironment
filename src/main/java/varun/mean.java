package varun;


import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Sampler;
import net.imglib2.algorithm.neighborhood.*;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import net.imglib2.algorithm.region.localneighborhood.*;
import net.imglib2.display.projector.sampler.IntervalSampler;
public class mean {
	
	
	
	public static < T extends RealType< T >> void mean( RandomAccessibleInterval<T> img, long radius){
		
		
		
	
		
	//	Random rand = new Random();
		// the number of dimensions
				long[] center = new long[img.numDimensions()];
				for (int d = 0; d < img.numDimensions(); ++d) {
					center[d] = (img.max(d) - img.min(d)) / 2 + img.min(d);
				}

				// define the center and radius
				Point centerdef = new Point(center);
				
				HyperSphere<T> testsphere = new HyperSphere<T>(img, centerdef, radius);

				Cursor<T> c = testsphere.cursor();

				long var = 0;

				while (c.hasNext()) {
					c.fwd();
					c.get().setReal(var);
					++var;

				}
				
				 Interval interval = Intervals.expand( img, -500 );
				
				 img=Views.interval(img, interval);
				 final Cursor < T > bound = Views.iterable(img).cursor();
				// long smvar=0;
				 
				 
				HyperSphere<T> smallsphere = new HyperSphere<T>(img, bound, radius/2);
				 Cursor<T> smallc = smallsphere.cursor();
				// float intensity = rand.nextFloat();
				 while(smallc.hasNext()){
					 
					 smallc.fwd();
					 
					 smallc.get().set(bound.get());
					 
					 
	//				

		//				while (smallc.hasNext()) {

			//				smallc.fwd();
							
				//			if (smvar > smallc.get().getRealDouble()) {

//								T value = smallc.get(); //Gives the value at that pixel
					//			smallc.get().setReal(smvar);
			//					value.setReal(smvar); // Sets the value

							

						//		smallc.get().setReal(smvar);
								
							//	++smvar;
						//	}

						//	}
				 }
				 
				//RandomAccessibleInterval<T> viewSource = Views.offsetInterval( img, new long[] { 0, 0 }, new long[]{ 100, 100 } );
				
				
		
		
	}
	
	
	
	
	public static void main (String[] args){
		
		
		final Img <FloatType> img= ArrayImgs.floats(1000,1000);
		
		long radi=400;
		mean(img, radi);
		ImageJFunctions.show(img);
		
		
	}

}
