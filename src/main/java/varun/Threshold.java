package varun;



import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;

import net.imglib2.view.Views;


public class Threshold {

	
	
	public static void createBitimage (RandomAccessibleInterval< FloatType > img, RandomAccessibleInterval<BitType > imgout, FloatType ThresholdValue){
		
		final Cursor< FloatType > bound =Views.iterable(img).localizingCursor();
		
		final RandomAccess< BitType > outbound = imgout.randomAccess();
		
		while(bound.hasNext()){
			
			bound.fwd();
			
			outbound.setPosition(bound);
			
			
			if(bound.get().compareTo(ThresholdValue) > 0){
				
				
				
				outbound.get().setOne();
				
			}
			
			
			else {
				
				outbound.get().setZero();
				
			}
		
			
			
			
		}
		
		
	}
	
	
	public static void main(String[] args){
		
		
		final Img< FloatType > img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		final Img< BitType>imgout= new ArrayImgFactory<BitType>().create(img, new BitType());
		
		 FloatType val= new FloatType(100);
		
		
		createBitimage(img,imgout,val);
		
		ImageJFunctions.show(imgout);
		
		
		
	}
	
	
	
	
	
	
}
