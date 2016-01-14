package varun;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class distanceTransform {

	public static void computeDistance(RandomAccessibleInterval<BitType> img,
			RandomAccessibleInterval<BitType> imgout) {

		final Cursor<BitType> bound = Views.iterable(img).cursor();

		final RandomAccess<BitType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			
			
			if (bound.get().getInteger() == 0) {

				double distance=0;
				double mindistance=0;
				double Huge=10;
				
				
				final Cursor<BitType> second = Views.iterable(img).cursor();

				while (second.hasNext()) {

					second.fwd();

					if (second.get().getInteger() == 1) {
						
						
						
						for(int i=0; i<second.numDimensions();++i){
							
						distance+=	Math.abs(bound.getIntPosition(i)-second.getIntPosition(i));
						}
						
					//	System.out.println(distance);

						// Compute mindistance somehow!

					}
					
					mindistance=Math.min(distance, Huge);

				}

//				mindistance=??
				
				outbound.get().setReal(mindistance);
				
			}

			
			
		}

	}

	public static void main(String[] args) {
		
		//final Img< FloatType > img=ArrayImgs.floats(100,100);
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		
		final Img< BitType > bitimg = Thresholder.threshold( img, new FloatType( 20 ), true, 1 );
		
		final Img< BitType>bitimgout= new ArrayImgFactory<BitType>().create(img, new BitType());
		
		
		
		
		computeDistance(bitimg,bitimgout);
		ImageJFunctions.show(img);
		//ImageJFunctions.show(bitimg);
		ImageJFunctions.show(bitimgout);
		

	}

}
